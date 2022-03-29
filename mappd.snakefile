
"""
MAPPD: Metagenomic Analysis Pipeline for Pathogen Discovery
Pipeline for finding pathogen DNA/RNA in metagenomic data
Includes trimming, assembly and annotation
Different pipelines can be selected and run in parallel
by altering the config.yaml file
"""

import os, sys

configfile: "config.yaml"

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
