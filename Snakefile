
"""
MAPPD: Metagenomic Analysis Pipeline for Pathogen Discovery
Pipeline for finding pathogen DNA/RNA in metagenomic data
Includes trimming, assembly and annotation
"""

configfile: "config.yaml"

include: "rules/preprocessing.rules"
include: "rules/assembly.rules"

rule all:
    input:
        expand("03_spades/{sample}", sample=config["samples"])