"""
Consume all benchmark files in the benchmark directory
Produce figures of time and memory usage for each rule
NOTE: had to put a dummy file as input to 'summarise_benchmarks
This ensures that benchmarks are calculated at the end of the pipeline
"""

rule summarise_benchmarks:
    input:
        # need to 'trick' snakemake into only running this at the end
        # but before the report is generated
        # ideally this will be, for example,  "annotation_summary.txt"
        expand("{sample}_host_depleted_1P.fastq", sample=config["samples"])
    output:
        "benchmarks/benchmarks.summary"
    shell:
        """
        {config[program_dir]}/scripts/summarise_benchmarks.py \
            -b benchmarks/ \
            -o {output}
        """

rule plot_bench_time:
    input:
        "benchmarks/benchmarks.summary"
    output:
        pdf = "benchmarks/bench_time.pdf",
        png = "benchmarks/bench_time.png"
    shell:
        """
        Rscript {config[program_dir]}/scripts/plot_benchmarks_time.R \
        {input} {output.pdf} {output.png}
        """

rule plot_bench_mem:
    input:
        "benchmarks/benchmarks.summary"
    output:
        pdf = "benchmarks/bench_mem.pdf",
        png = "benchmarks/bench_mem.png"
    shell:
        """
        Rscript {config[program_dir]}/scripts/plot_benchmarks_mem.R \
        {input} {output.pdf} {output.png}
        """

rule draw_dag:
    input:
        "mappd.snakefile"
    output:
        "benchmarks/dag.png"
    shell:
        "snakemake -s {input} --rulegraph 2> /dev/null | dot -T png > {output}"
