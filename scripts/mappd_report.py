#!/usr/bin/env python

"""
Generate a html report for the mappd pipeline
Will generate a report specific to the files provided
For example, if the pipeline is run with only the trinity assembler,
this report will automatically only include results from
that assembly
"""

import sys, os
from snakemake.report import data_uri_from_file
import maketable

# this way the report should work no matter what pipeline components are run
# and I'll also be able to refer to them easily

def generate_report(config_file="", dag_graph="",
                    bench_time="", bench_mem="",
                    trim_summary="",
                    host_table="",
                    mapping_figure="",
                    spades_assembly="", spades_bandage="",
                    trinity_assembly="", trinity_bandage="",
                    spades_diamond="",
                    trinity_diamond="",
                    spades_blastn="",
                    trinity_blastn=""):

    report = """
==========================================================================================
MAPPD: Metagenomic Analysis Pipeline for Pathogen Discovery
==========================================================================================

"""
    if dag_graph:
        report += """
DAG graph of pipeline components
================================
    A Directed Acyclic Graph (DAG) graph of the steps carried out in this pipeline is given below

"""
        report += "\t.. image:: " + data_uri_from_file(dag_graph)[0] + "\n"

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
        report += "\t.. image:: " + data_uri_from_file(bench_mem)[0] + "\n"

    if trim_summary:
        report += """
Trimming
=================
    The reads were trimmed using Trimmomatic.
    The figure below shows how many read pairs were both retained, how many single pairs were created
    and how many reads were dropped altogether.

"""
        report += "\t.. image:: " + data_uri_from_file(trim_summary)[0] + "\n"

    if host_table:
        host_string = maketable.make_table_from_csv(host_table, sep="\t")
        # also want to extract out the text for the host species downloaded
        # will read the top 10 table into memory and extract those lines
        host_lines = [line for line in open(host_table)]
        # this ignores the header, extracts number of lines / hosts according to config
        host_strings = host_lines[1:config_file["hosts_to_download"]+1]
        # this gets just the species name from the host string
        hosts = [host.split("\t")[3] for host in host_strings]

        report += """
Ribosomal RNA and host removal
===============================
*rRNA mapping*
------------------
    The trimmed reads were aligned to the SILVA ribosomal RNA databases, including the Long Sub Unit (LSU)
    and Small Sub Unit (SSU) categories, and matching reads were removed.

*Host mapping*
---------------
    To identify the most abundant organism (potentially the 'host'),
    a subset of 200,000 putative mRNA sequences from each sample were extracted and
    individually assembled using SPAdes.
    The top 10 contigs larger than 1,000 bps from the assemblies were blasted
    against the NCBI nr database and the most abundant, best matches were used
    to assign the most likely host species as {}.
    All genetic data from was then extracted from the NCBI nr database
    and was used to identify host sequences in the samples.
    
""".format(hosts[0])
        report += host_string + "\n"
        report += "\t.. image:: " + data_uri_from_file(mapping_figure)[0] + "\n"


    if spades_bandage:
        report += """
Assembly
=================
    The reads were assembled using SPAdes.
    The figure below gives a representation of the scaffolds with at least 10x coverage

"""
        # NOTE: spades_bandage is a 'named list' due to wildcard expansion
        # thus, have to take first element of list to get str for data_uri
        report += "\t.. image:: " + data_uri_from_file(spades_bandage[0])[0] + "\n"

    if trinity_bandage:
        report += """
    The reads were assembled using Trinity.
    The figure below gives a representation of the scaffolds with at least 10x coverage

"""
        # NOTE: spades_bandage is a 'named list' due to wildcard expansion
        # thus, have to take first element of list to get str for data_uri
        report += "\t.. image:: " + data_uri_from_file(trinity_bandage[0])[0] + "\n"

    return(report)