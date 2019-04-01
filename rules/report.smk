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

sys.path.insert(0, config["program_dir"] + "/scripts") # this allows me to import modules in this folder

from snakemake.utils import report
from mappd_report import generate_report

rule full_run_report:
    input:
        dag_graph = "benchmarks/dag.png",
        bench_time = "benchmarks/bench_time.png",
        bench_mem = "benchmarks/bench_mem.png",
        trim_summary = "logs/" + config["sub_dirs"]["trim_dir"] + "/trim_summary.png",
        # NOTE: the below parameters are received as a 'named list' due to wildcard expansion
        spades_assembly = expand(config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts.fasta",
        sample=config["samples"]),
        spades_bandage = expand(config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/assembly_graph_10x.png",
            sample=config["samples"]),
        trinity_bandage = expand(config["sub_dirs"]["assembly_dir"] + "/trinity/{sample}_trinity/{" \
            "sample}_trinity.Trinity.png", sample=config["samples"])

    output:
        "full_report.html"
    run:
        sphinx_str = generate_report(config_file=config, dag_graph=input.dag_graph,
                                     bench_mem=input.bench_mem, bench_time=input.bench_time,
                                     trim_summary=input.trim_summary,
                                     spades_bandage=input.spades_bandage,
                                     trinity_bandage=input.trinity_bandage)
        report(sphinx_str, output[0], metadata="Author: Matthew Neave (matthew.neave@csiro.au)")
