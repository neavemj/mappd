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

rule plot_bench_time:
    input:
        "benchmarks/benchmarks.summary"
    output:
        "benchmarks/bench_time.pdf"
    shell:
         """
         Rscript {config[program_dir]}/scripts/plot_benchmarks_time.R \
         {input} {output} benchmarks/bench_time.png
         """

rule plot_bench_mem:
    input:
         "benchmarks/benchmarks.summary"
    output:
          "benchmarks/bench_mem.pdf"
    shell:
         """
         Rscript {config[program_dir]}/scripts/plot_benchmarks_mem.R \
         {input} {output} benchmarks/bench_mem.png
         """
