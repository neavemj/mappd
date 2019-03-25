"""
Consume all benchmark files in the benchmark directory
Produce figures of time and memory usage for each rule
"""

rule summarise_benchmarks:
    input:
        "benchmarks/"
    output:
        "benchmarks/benchmarks.summary"
    shell:
        """
        scripts/summarise_benchmarks.py \
            -b {input} \
            -o {output}
        """

rule plot_benchmarks:
    input:
        "benchmarks/benchmarks.summary"
    output:
        "benchmarks/bench_summary.pdf"
    shell:
         """
         Rscript {config[program_dir]}/scripts/plot_benchmarks.R \
         {input} {output} benchmarks/bench_summary.png
         """
