
"""
MAPPD: Metagenomic Analysis Pipeline for Pathogen Discovery
Pipeline for finding pathogen DNA/RNA in metagenomic data
Includes trimming, assembly and annotation
"""

configfile: "config.yaml"

include: "rules/trim_quality.rules"

rule all:
    input:
        "trim/simulated_sample.1m.trimmed.fastq"