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
    input:
        getFastq
    output:
        R1_P = "01_trimmomatic/{sample}_1P.fastq.gz",
        R1_U = "01_trimmomatic/{sample}_1U.fastq.gz",
        R2_P = "01_trimmomatic/{sample}_2P.fastq.gz",
        R2_U = "01_trimmomatic/{sample}_2U.fastq.gz"
    params:
        qual = config["trimmomatic_quality"],
        adapters = config["program_dir"] + config["trimmomatic_adapters"],
        minlen = config["trimmomatic_minlen"]
    threads: 8
    log:
        "logs/01_trimmomatic/{sample}.log"
    benchmark:
        "benchmarks/01_trimmomatic/{sample}.txt"
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
        expand("logs/01_trimmomatic/{sample}.log", sample=config["samples"])
    output:
        "logs/01_trimmomatic/logs.summary"
    shell:
        """
        {config[program_dir]}/scripts/summarise_trimmomatic.py \
        -i {input} -o {output}
        """

rule plot_trimmomatic_results:
    input:
        "logs/01_trimmomatic/logs.summary"
    output:
        "logs/01_trimmomatic/trim_summary.pdf"
    shell:
        """
        Rscript {config[program_dir]}/scripts/plot_trim.R \
        {input} {output} logs/01_trimmomatic/trim_summary.png
        """
