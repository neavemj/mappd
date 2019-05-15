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

def generate_report(config="", dag_graph="",
                    bench_time="", bench_mem="",
                    software_versions="",
                    overall_figure = "",
                    spades_assembly="", spades_bandage="",
                    trinity_assembly="", trinity_bandage="",
                    spades_diamond="",
                    trinity_diamond="",
                    spades_blastn="",
                    trinity_blastn="",
                    ):

    # if a superkingdom is not detected, the figs and tables won't be created
    # e.g. if no viruses in a sample, there won't be any figures for these
    # will go through and check which ones exist
    potential_taxa_files = {
        "euk_figure": config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance_top10.euk.png",
        "euk_table": config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance_top10.euk.tsv",
        "bac_figure": config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance_top10.bac.png",
        "bac_table": config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance_top10.bac.tsv",
        "vir_figure": config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance_top10.vir.png",
        "vir_table": config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance_top10.vir.tsv",
    }

    taxa_files = {}
    for taxa in potential_taxa_files:
        if os.path.isfile(potential_taxa_files[taxa]):
            taxa_files[taxa] = potential_taxa_files[taxa]

    report = """

.. contents:: MAPPD: Metagenomic Analysis Pipeline for Pathogen Discovery

|
|

_______

    """

    report += """

1   Introduction
==================
MAPPD is a general pipeline for the identification of organisms in a metagenomic sample,
although it is targeted toward the identification of pathogens.
The pipeline uses a strategy of read quality trimming, host identification,
read assembly, and annotation using various blast and diamond searches.
MAPPD does not require prior information about the sample (e.g. host),
as this information is determined by classifying read sub-sets.

**Important**
---------------
Metagenomic analysis can be a useful technology for screening samples in cases
where a pathogen is unknown. However, the classication of sequence fragments
based on the highest identity in a database does not necessarily mean that a
pathogen is present, only that this is the 'best' match. This report provides
the percent identity of database hits and the location of the particular contigs.
It may be necessary to check important classifications manually.
Additional lab-based tests are required to confirm pathogen identification.

|
|

_________

2   Technical Summary
========================
The software verions used in this pipeline are given below.

"""
    software_string = maketable.make_table_from_csv(software_versions, sep="\t")
    report += software_string + "\n"

    report += """

|
|

_________

3   Data Quality and Overall Classifications
===============================================
The raw data were trimmed for quality and adapters using `Trimmomatic`_.
The cleaned reads were then aligned to the SILVA ribosomal RNA databases,
including the Long Sub Unit (LSU) and Small Sub Unit (SSU) categories,
and matching reads were removed.

The remaining reads were then classified using iterative assemblies,
blasts and diamond searches. Reads that could not be classified after these
processes are shown as the grey 'Unannotated' bar below.

.. _Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic

    """
    report += "\t.. image:: " + data_uri_from_file(overall_figure)[0] + "\n"

    if spades_bandage:
        report += """

|
|

________

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

    if "euk_figure" in taxa_files:
        report += """

|
|

________

4   Summary Classifications
=============================
*Eukaryotes*
---------------
The contigs were annotated using Diamond blastx.
The figure most abundant Eukaryotic organisms
in all of the samples combined.

"""
        report += "\t.. image:: " + data_uri_from_file(taxa_files["euk_figure"])[0] + "\n"

        report += """

|

The table shows how many reads were assigned to each `organism`_.

.. _organism: 04_annotation/diamond/diamond_blastx_abundance.euk

"""

        euk_string = maketable.make_table_from_csv(taxa_files["euk_table"], sep="\t")
        report += euk_string + "\n"

    if "bac_figure" in taxa_files:
        report += """

|

*Bacteria*
---------------
The figure below shows the most abundant Bacterial organisms
in all of the samples combined.

"""
        report += "\t.. image:: " + data_uri_from_file(taxa_files["bac_figure"])[0] + "\n"

        report += """

|

The table shows how many reads were assigned to each bacteria.

"""

        bac_string = maketable.make_table_from_csv(taxa_files["bac_table"], sep="\t")
        report += bac_string + "\n"

    if "vir_figure" in taxa_files:
        report += """

|

*Viruses*
---------------
The figure below shows the most abundant Viral organisms
in all of the samples combined.

"""
        report += "\t.. image:: " + data_uri_from_file(taxa_files["vir_figure"])[0] + "\n"

        report += """

|

The table shows how many reads were assigned to each virus.

"""

        vir_string = maketable.make_table_from_csv(taxa_files["vir_table"], sep="\t")
        report += vir_string + "\n"

    if dag_graph:
        report += """

|
|

__________

Detailed DAG graph of pipeline structure
===========================================
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

    return(report)