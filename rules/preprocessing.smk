"""
These rules will trim sequencing reads

The inputs are:
    - raw Illumina sequencing reads
The outputs are:
    - cleaned Illumina sequencing reads
"""

# for some reason I have put this function here to get wildcards to work
# it wont work if directly in the rule
def getFastq(wildcards):
    return config['samples'][wildcards.sample]


# TODO: might need a trimmomatic SE mode
rule trimmomatic_PE:
    message:
        """
        Trimming for quality and Illumina adapters using Trimmomatic
        """
    input:
        getFastq
    output:
        R1_P = config["sub_dirs"]["trim_dir"] + "/{sample}_1P.fastq.gz",
        R1_U = config["sub_dirs"]["trim_dir"] + "/{sample}_1U.fastq.gz",
        R2_P = config["sub_dirs"]["trim_dir"] + "/{sample}_2P.fastq.gz",
        R2_U = config["sub_dirs"]["trim_dir"] + "/{sample}_2U.fastq.gz"
    params:
        qual = config["trimmomatic_quality"],
        adapters = config["program_dir"] + config["trimmomatic_adapters"],
        minlen = config["trimmomatic_minlen"]
    threads: 8
    log:
        "logs/trimmomatic_PE/{sample}.log"
    benchmark:
        "benchmarks/trimmomatic_PE/{sample}.txt"
    shell:
        """
        trimmomatic PE \
            -threads {threads} \
            {input} {output.R1_P} {output.R1_U} {output.R2_P} {output.R2_U} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:{params.qual} MINLEN:{params.minlen} \
            2> {log}
        """

rule summarise_trimmomatic_log:
    input:
        expand("logs/" + config["sub_dirs"]["trim_dir"] + "/{sample}.log",
            sample=config["samples"])
    output:
        "logs/" + config["sub_dirs"]["trim_dir"] + "/logs.summary"
    shell:
        """
        {config[program_dir]}/scripts/summarise_trimmomatic.py \
        -i {input} -o {output}
        """

rule plot_trimmomatic_results:
    input:
        "logs/" + config["sub_dirs"]["trim_dir"] + "/logs.summary"
    output:
        pdf = "logs/" + config["sub_dirs"]["trim_dir"] + "/trim_summary.pdf",
        png = "logs/" + config["sub_dirs"]["trim_dir"] + "/trim_summary.png"
    shell:
        """
        Rscript {config[program_dir]}/scripts/plot_trim.R \
        {input} {output.pdf} {output.png}
        """
