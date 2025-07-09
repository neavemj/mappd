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
        config["sub_dirs"]["annotation_dir"] + "/diamond/png_file_names.txt",
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
        """
		module load graphviz/2.47.0
		
		snakemake -s {input} --rulegraph 2> /dev/null | dot -T png > {output}
		"""

rule record_end:
    input:
        # need to 'trick' snakemake into only running this at the end
        # ideally this will be, for example,  "annotation_summary.txt"
        start = "logs/start_time.txt",
        dummy = config["sub_dirs"]["annotation_dir"] + "/diamond/png_file_names.txt",
    output:
        end = "logs/end_time.txt",
        times = "logs/time_log.txt"
    shell:
        """
        echo -e "end time\t"$(date) > {output.end}
        cat {input.start} {output.end} > {output.times}
        """

rule get_package_versions:
    input:
        config["program_dir"] + "config/software_list.txt"
    output:
        "logs/software_versions.txt"
    params:
        # for some reason, I have to define the sed pattern here, then pass it in
        # otherwise it is weirdly expanded in the actual shell command
        sed_pat1 = r"s/\(.*\)/^\1\t/g",
        sed_pat2 = r"s/ \+/\t/g"
    shell:
        # lists all packages in the conda environment
        # replaces large whitespace with tabs
        # greps for specific packages with ^ and \t to ensure complete match
        # cuts just the first 2 columns of interest
        """
        conda list | \
            sed "{params.sed_pat2}" | \
            grep -f <(cat {input} | sed "{params.sed_pat1}") | \
            cut -f 1,2 > {output}
        """

rule compile_technical_summary:
    input:
        start = "logs/start_time.txt",
        end = "logs/end_time.txt",
        software = "logs/software_versions.txt",
        config = rules.create_sample_config.output[0]
    output:
        "logs/technical_summary.ReST"
    shell:
        # note: will pass sample info directly in from the config file
        """
        {config[program_dir]}/scripts/technical_ReST.py \
            -s {input.start} \
            -e {input.end} \
			-c {input.config} \
            -f {input.software} \
            -o {output}
        """






