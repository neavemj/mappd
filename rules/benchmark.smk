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
        "snakemake -s {input} --rulegraph 2> /dev/null | dot -T png > {output}"

rule record_end:
    input:
        # need to 'trick' snakemake into only running this at the end
        # ideally this will be, for example,  "annotation_summary.txt"
        config["sub_dirs"]["annotation_dir"] + "/diamond/png_file_names.txt",
    output:
        "logs/end_time.txt"
    shell:
        """
        echo "End time\t"$(date) > {output}
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
        # Adds 'version' in front of the software
        """
        conda list | \
            sed "{params.sed_pat2}" | \
            grep -f <(cat {input} | sed "{params.sed_pat1}") | \
            cut -f 1,2 | \
            sed "s/^/Version /g" >> {output}
        """

rule get_sample_list:
    output:
        "logs/sample_list.txt"
    params:
        # had to put all the patterns here to ensure
        # snakemake doesn't mangle them
        # see below for explanation of patterns
        sed1 = r"s/], /\n/g",
        sed2 = r"s/{\|}\|\[\|\]//g",
        sed3 = r"s/^/Sample /g",
        sed4 = r"s/:/\t/g",
    shell:
        # prints a dictionary of samples in the config file
        # creates a newline for each sample
        # removes square brackets and braces
        # adds 'Sample' to the beginning of each line
        # replaces the colon with a tab for the technical summary table
        """
        echo {config[samples]} | \
            sed "{params.sed1}" | \
            sed "{params.sed2}" | \
            sed "{params.sed3}" | \
            sed "{params.sed4}" > {output}
        """

rule compile_technical_summary:
    input:
        "logs/start_time.txt",
        "logs/end_time.txt",
        "logs/sample_list.txt",
        "logs/software_versions.txt",
    output:
        "logs/technical_summary.txt"
    shell:
        """
        echo "Parameter\tValue" > {output} &&
        cat {input} >> {output}
        """










