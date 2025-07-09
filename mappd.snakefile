
"""
MAPPD: Metagenomic Analysis Pipeline for Pathogen Discovery
Pipeline for finding pathogen DNA/RNA in metagenomic data
Includes trimming, assembly and annotation
Different pipelines can be selected and run in parallel
by altering the config.yaml file
"""

import os, sys
from os.path import join
from snakemake.io import glob_wildcards
configfile: "config.yaml"

sys.path.insert(0, config["program_dir"] + "scripts")
import shared_vars


# getting the list of samples here so it's available to all rules

RAW_DIR = config["raw_dir"]

SAMPLES = shared_vars.get_samples(RAW_DIR, config["illumina_machine"])

rules_dir = os.path.join(os.path.expanduser(config["program_dir"]), "rules")

include: os.path.join(rules_dir, "preprocessing.smk")
include: os.path.join(rules_dir, "rRNA_depletion.smk")
include: os.path.join(rules_dir, "host_depletion.smk")
include: os.path.join(rules_dir, "assembly.smk")
include: os.path.join(rules_dir, "annotation.smk")
include: os.path.join(rules_dir, "benchmark.smk")
include: os.path.join(rules_dir, "report.smk")


rule all:
    input:
        #expand("{pipe}_report.html", pipe=config["pipeline"])
        "mappd_report.html"
