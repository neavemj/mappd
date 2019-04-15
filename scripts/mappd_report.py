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

def generate_report(config_file="config", dag_graph="",
                    bench_time="", bench_mem="",
                    trim_summary="",
                    LSU_table="",
                    LSU_figure="",
                    SSU_table="",
                    SSU_figure="",
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

    if LSU_table:
        report += """
Ribosomal RNA removal
======================
    The trimmed reads were aligned to the SILVA databases, including the Long Sub Unit (LSU)
    and Small Sub Unit (SSU) categories.
    Matching reads were used to identify the host species for the host depletion step
    and were then removed from further analysis.

*LSU mapping*
---------------
    The top 10 organisms with the most hits to the LSU database are summarised below

"""
        LSU_string = maketable.make_table_from_csv(LSU_table[0], sep="\t").replace(";", " ")
        report += LSU_string + "\n"
        report += "\t.. image:: " + data_uri_from_file(LSU_figure)[0] + "\n"

        report += """
*SSU mapping*
---------------
    The top 10 organisms with the most hits to the SSU database are summarised below

"""
        SSU_string = maketable.make_table_from_csv(SSU_table[0], sep="\t").replace(";", " ")
        report += SSU_string + "\n"
        report += "\t.. image:: " + data_uri_from_file(SSU_figure)[0] + "\n"

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