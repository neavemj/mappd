
"""
MAPPD: Metagenomic Analysis Pipeline for Pathogen Discovery
Pipeline for finding pathogen DNA/RNA in metagenomic data
Includes trimming, assembly and annotation
"""

import os

configfile: 'config.yaml'

rules_dir = os.path.join(os.path.expanduser(config['program_dir']), 'rules')

include: os.path.join(rules_dir, 'preprocessing.smk')
include: os.path.join(rules_dir, 'assembly.smk')
include: os.path.join(rules_dir, 'benchmark.smk')

rule all:
    input:
        #expand("01_trimmomatic/{sample}_1P.fastq.gz", sample=config['samples'])
        #expand("02_spades/{sample}", sample=config["samples"])
        "benchmarks/bench_summary.pdf"