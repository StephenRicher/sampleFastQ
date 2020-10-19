#!/usr/bin/env python3

import os
import re
import sys
import tempfile
import pandas as pd
from snake_setup import set_config, load_samples

BASE = workflow.basedir

# Define path to conda environment specifications
ENVS = f'{BASE}/workflow/envs'
# Defne path to custom scripts directory
SCRIPTS = f'{BASE}/workflow/scripts'

if not config:
    configfile: f'{BASE}/config/config.yaml'

# Defaults configuration file - use empty string to represent no default value.
default_config = {
    'workdir':           workflow.basedir,
    'tmpdir':            tempfile.gettempdir(),
    'threads':           workflow.cores,
    'data':              ''          ,
    'paired':            ''          ,
    'spliced':           ''          ,
    'genome':
        {'build':          'sequence',
         'regions':        ''        ,
         'index':          None      ,
         'sequence':       ''        ,},
    'cutadapt':
        {'forwardAdapter': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
         'reverseAdapter': 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
         'overlap':         3                                 ,
         'errorRate':       0.1                               ,
         'minimumLength':   0                                 ,
         'qualityCutoff':  '0,0'                              ,
         'GCcontent':       50                                ,},
}

config = set_config(config, default_config)

workdir: config['workdir']

# Read path to samples in pandas
samples = load_samples(config['data'], config['paired'])

wildcard_constraints:
    single = '|'.join(samples['single']),
    sample = '|'.join(samples['sample']),
    rep = '|'.join(samples['rep'].unique()),

rule all:
    input:
        expand('fastq/{single}.fastq.gz', single=list(samples["single"]))


def cutadaptOutput():
    if config['paired']:
        return ['dat/fastq/trimmed/{sample}-R1.trim.fastq.gz',
                'dat/fastq/trimmed/{sample}-R2.trim.fastq.gz']
    else:
        return ['fastq/trimmed/{sample}.trim.fastq.gz']


def cutadaptCmd():

    if config['paired']:
        cmd = ('cutadapt -a {params.forwardAdapter} -A {params.reverseAdapter} '
               '-o {output.trimmed[0]} -p {output.trimmed[1]} ')
    else:
        cmd = 'cutadapt -a {params.forwardAdapter} -o {output.trimmed[0]} '
    cmd += ('--overlap {params.overlap} --error-rate {params.errorRate} '
            '--minimum-length {params.minimumLength} '
            '--quality-cutoff {params.qualityCutoff} '
            '--gc-content {params.GCcontent} '
            '--cores {threads} {input} > {output.qc} 2> {log}')
    return cmd


rule cutadapt:
    input:
        lambda wc: samples.xs(wc.sample, level=1)['path']
    output:
        trimmed = cutadaptOutput(),
        qc = 'qc/cutadapt/unmod/{sample}.cutadapt.txt'
    params:
        forwardAdapter = config['cutadapt']['forwardAdapter'],
        reverseAdapter = config['cutadapt']['reverseAdapter'],
        overlap = config['cutadapt']['overlap'],
        errorRate = config['cutadapt']['errorRate'],
        minimumLength = config['cutadapt']['minimumLength'],
        qualityCutoff = config['cutadapt']['qualityCutoff'],
        GCcontent = config['cutadapt']['GCcontent']
    log:
        'logs/cutadapt/{sample}.log'
    conda:
        f'{ENVS}/cutadapt.yaml'
    threads:
        config["threads"]
    shell:
        cutadaptCmd()


if config['spliced']:

    def hisat2Cmd():
        if config['paired']:
            return ('hisat2 -x {params.index} -1 {input.reads[0]} '
                '-2 {input.reads[1]} --threads {threads} '
                '--summary-file {output.qc} > {output.sam} 2> {log}')
        else:
            return ('hisat2 -x {params.index} -p {threads} '
                '-U {input.reads[0]} --summary-file {output.qc} '
                '> {output.sam} 2> {log}')

    rule hisat2:
        input:
            reads = rules.cutadapt.output.trimmed
        output:
            sam = pipe('mapped/{sample}.sam'),
            qc = 'qc/{sample}.hisat2.txt'
        params:
            index = config['genome']['index']
        group:
            'map'
        log:
            'logs/hisat2/{sample}.log'
        conda:
            f'{ENVS}/hisat2.yaml'
        threads:
            max(1, int(config['threads'] / 2))
        shell:
            hisat2Cmd()

else:

    rule bowtie2Build:
        input:
            config['genome']['sequence']
        output:
            expand('dat/genome/index/{build}.{n}.bt2',
                   n=['1', '2', '3', '4', 'rev.1', 'rev.2'],
                   build=config["genome"]["build"])
        params:
            basename = f'dat/genome/index/{config["genome"]["build"]}'
        log:
            'logs/bowtie2Build.log'
        conda:
            f'{ENVS}/bowtie2.yaml'
        threads:
            config["threads"]
        shell:
            'bowtie2-build --threads {threads} {input} {params.basename} &> {log}'


    def bowtie2Cmd():
        if config['paired']:
            return (
                'bowtie2 -x {params.index} -1 {input.reads[0]} '
                '-2 {input.reads[1]} --threads {threads} '
                '> {output.sam} 2> {log}; cp {log} {output.qc}')
        else:
            return (
                'bowtie2 -x {params.index} -U {input.reads[0]} '
                '--threads {threads} '
                '> {output.sam} 2> {log}; cp {log} {output.qc}')

    rule bowtie2:
        input:
            reads = rules.cutadapt.output.trimmed,
            index = rules.bowtie2Build.output
        output:
            sam = pipe('mapped/{sample}.sam'),
            qc = 'qc/bowtie2/{sample}.bowtie2.txt'
        params:
            index = f'dat/genome/index/{config["genome"]["build"]}'
        group:
            'map'
        log:
            'logs/bowtie2/{sample}.log'
        conda:
            f'{ENVS}/bowtie2.yaml'
        threads:
            max(1, int(config['threads'] / 2))
        shell:
            bowtie2Cmd()


def getMapOutput(wc):
    if config['spliced']:
        return rules.hisat2.output.sam
    else:
        return rules.bowtie2.output.sam


rule sortBAM:
    input:
        getMapOutput
    output:
        'mapped/{sample}.sort.bam'
    group:
        'map'
    log:
        'logs/sortBAM/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        max(1, int(config['threads'] / 2))
    shell:
        'samtools sort -@ {threads} {input} > {output} 2> {log}'


rule indexBAM:
    input:
        rules.sortBAM.output
    output:
        f'{rules.sortBAM.output}.bai'
    log:
        'logs/indexBAM/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        config['threads']
    shell:
        'samtools index -@ {threads} {input} &> {log}'


rule filterRegions:
    input:
        bam = rules.sortBAM.output,
        index = rules.indexBAM.output
    output:
        'mapped/{sample}.filt.bam'
    params:
        bed = config['genome']['regions']
    log:
        'logs/filterRegions/{sample}.log'
    conda:
        f'{ENVS}/samtools.yaml'
    threads:
        max(1, config['threads'])
    shell:
        'samtools view -b -@ {threads} -L {params.bed} '
        '{input.bam} > {output} 2> {log}'


rule getFastQIDs:
    input:
        rules.filterRegions.output.bam
    output:
        'fastq/{sample}-validIDs.txt'
    log:
        'logs/getFastQIDs/{sample}.log'
    threads:
        max(1, config['threads'])
    conda:
        f'{ENVS}/samtools.yaml'
    shell:
        '(samtools view -@ {threads} {input} | awk -f {SCRIPTS}/getID.awk '
        '| sed s"/\/[12]$//" | uniq > {output}) 2> {log}'


rule filterFastQ:
    input:
        ids = rules.getFastQIDs.output,
        reads = lambda wc: samples.xs(f'{wc.sample}-{wc.read}', level=2)['path']
    output:
        'fastq/{sample}-{read}.fastq.gz'
    log:
        'logs/filterFastQ/{sample}-{read}.log'
    shell:
        'grep -A 3 -Ff {input.ids} <(zcat -f {input.reads}) '
        '| grep -v "^--$" | gzip > {output} 2> {log} '
