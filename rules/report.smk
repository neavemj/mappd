"""
mappd reporting rule

This rule defines the pipeline components by the input files that it requires
I.E., snakemake knows which rules to run by the input files that are required here
The generate_report function contains named parameters
so that 'if' statements can be used to report only items run in the pipeline
The mappd_report.py script should NOT need modification for new pipelines
To create a new pipeline, create a new rule below and include the steps you want
in the input files, then give only these files to the generate_report function
The rule for each pipeline should end with '_report.html',
then the pipeline name can be specified in the config.yaml,
and will be included as a wildcard in 'mappd.snakefile'
This allows multiple pipelines to be run at once.
"""

import sys
sys.path.insert(0, config["program_dir"] + "/scripts") # this allows me to import modules in this folder
from snakemake.utils import report
from mappd_report import generate_report


rule test_report:
    input:
        dag_graph = "benchmarks/dag.png",
        #bench_time = "benchmarks/bench_time.png",
        #bench_mem = "benchmarks/bench_mem.png",
        software_versions = "logs/software_versions.txt",
        trim_summary = "logs/trimmomatic_PE/trim_summary.png",
        host_table = config["sub_dirs"]["depletion_dir"] + "/host/largest_contigs.blastn.tax.wide",
        mapping_figure = "logs/mapping_summary.png",
        euk_figure = config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance.euk.png",
        bac_figure = config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance.bac.png",
        vir_figure = config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance.vir.png",
        report_css = config["program_dir"] + "config/report.css"
    output:
        "test_report.html"
    run:
        sphinx_str = generate_report(config_file=config, dag_graph=input.dag_graph,
                                     #bench_mem=input.bench_mem, bench_time=input.bench_time,
                                     software_versions=input.software_versions,
                                     trim_summary=input.trim_summary,
                                     host_table=input.host_table,
                                     mapping_figure=input.mapping_figure,
                                     euk_figure=input.euk_figure,
                                     bac_figure=input.bac_figure,
                                     vir_figure=input.vir_figure,
                                     )
        report(sphinx_str, output[0], stylesheet=input.report_css, metadata="Author: Matthew Neave (matthew.neave@csiro.au)")

