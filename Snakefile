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
    'data':              ''          ,
    'paired':            ''          ,
    'genome':
        {'build':          'sequence',
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
        workflow.cores
    shell:
        'cutadapt -a {params.forwardAdapter} '
        '--overlap {params.overlap} --error-rate {params.errorRate} '
        '--minimum-length {params.minimumLength} '
        '--quality-cutoff {params.qualityCutoff} '
        '--gc-content {params.GCcontent} '
        '--cores {threads} -o {output.trimmed[0]} {input} '
        '> {output.qc} 2> {log}'


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
        workflow.cores
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

rule bowtie2Map:
    input:
        reads = rules.cutadapt.output.trimmed,
        index = rules.bowtie2Build.output
    output:
        sam = pipe('mapped/{sample}.sam'),
        qc = 'qc/bowtie2/{sample}.bowtie2.txt'
    params:
        index = f'genome/index/{config["genome"]["build"]}'
    group:
        'map'
    log:
        'logs/bowtie2Map/{sample}.log'
    conda:
        f'{ENVS}/bowtie2.yaml'
    threads:
        workflow.cores - 1
    shell:
        bowtie2Cmd()


rule getFastQIDs:
    input:
        rules.bowtie2Map.output.sam
    output:
        'fastq/{sample}-validIDs.txt'
    params:
        chr = 11
    group:
        'map'
    log:
        'logs/getFastQIDs/{sample}.log'
    shell:
        """awk "\$3=={params.chr} {{print "@"\$1}}" | uniq > {output} """


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