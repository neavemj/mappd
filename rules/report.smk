"""
mappd reporting rule
"""

from snakemake.utils import report
from scripts.mappd_report import generate_report

rule full_run_report:
    input:
        trim_summary = "logs/01_trimmomatic/trim_summary.png",
        bench_time = "benchmarks/bench_time.png",
        bench_mem = "benchmarks/bench_mem.png",
        spades_assembly = expand(directory("02_spades/{sample}"), sample=config["samples"])
    output:
        "full_report.html"
    run:
        sphinx_str = generate_report(config_file=config, bench_mem=input.bench_mem, bench_time=input.bench_time,
                                     trim_summary=input.trim_summary)
        report(sphinx_str, output[0])

