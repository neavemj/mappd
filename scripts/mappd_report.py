#!/usr/bin/env python

"""
Generate a html report for the mappd pipeline
Will generate a report specific to the files provided
For example, if the pipeline is run with only the trinity assembler,
this report will automatically only include results from
that assembly
"""

import os
from snakemake.report import data_uri_from_file

# this way the report should work no matter what pipeline components are run
# and I'll also be able to refer to them easily

def generate_report(config_file="config", bench_time="", bench_mem="", trim_summary="",
                    spades_assembly="", spades_bandage="", trinity_assembly="", spades_diamond="", trinity_diamond="",
                    spades_blastn="", trinity_blastn=""):
    print(bench_mem)
    report = """
==========================================================================================
MAPPD: Metagenomic Analysis Pipeline for Pathogen Discovery
==========================================================================================
    
"""

    if bench_time:
        report += """
Benchmarks
=================
    The time taken for each process, and each sample, is given below
    
"""
        report += "\t.. image:: " + data_uri_from_file(bench_time)[0] + "\n"

    if bench_mem:
        report += """
    The maximum memory required for each process, and each sample, is given below
    
"""
        report += "\t.. image:: " + data_uri_from_file(bench_mem, "temp.png")[0] + "\n"

    if trim_summary:
        report += """
Trimming
=================
    The reads were trimmed using Trimmomatic.
    The figure below shows how many read pairs were both retained, how many single pairs were created 
    and how many reads were dropped altogether.
    
"""
        report += "\t.. image:: " + data_uri_from_file(trim_summary)[0] + "\n"

    if spades_bandage:
        report += """
Assembly
=================
    The reads were assembled using Trimmomatic.
    The figure below gives an overview of the contigs with at least 10x coverage

"""
        report += "\t.. image:: " + data_uri_from_file(spades_bandage)[0] + "\n"

    return(report)