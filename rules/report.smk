"""
mappd reporting rule
"""

from snakemake.utils import report
#from scripts.mappd_report import generate_report

rule full_run_report:
    input:
        "logs/01_trimmomatic/trim_summary.png",
        "benchmarks/bench_time.png",
        "benchmarks/bench_mem.png"
    output:
        "full_report.html"
    shell:
        """
        ls {input} > {output}
        """

rule trinity_diamond_report:
    input:
        "logs/01_trimmomatic/trim_summary.png",
        "benchmarks/bench_time.png",
        "benchmarks/bench_mem.png"
    output:
        "trinity_diamond_report.html"
    shell:
        """
        ls {input} > {output}
        """
#    run:
#        sphinx_str = get_sphinx_report(config, input)
#        report(sphinx_str, output)
