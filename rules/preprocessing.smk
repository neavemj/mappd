"""
These rules will trim sequencing reads

The inputs are:
    - raw Illumina sequencing reads
The outputs are:
    - cleaned Illumina sequencing reads
"""

# record start time of the run

rule record_start:
    output:
        "logs/start_time.txt"
    shell:
        """
        echo -e "start time\t"$(date) > {output}
        """

# for some reason I have put this function here to get wildcards to work
# it wont work if directly in the rule
def getFastq(wildcards):
    return config['samples'][wildcards.sample]


# TODO: might need a trimmomatic SE mode
rule trimmomatic_PE:
    message:
        """
        ** preprocessing **
        Trimming {wildcards.sample} for quality and Illumina adapters using Trimmomatic
        """
    input:
        reads = getFastq,
        # trick to get date recorded at this first step
        date = "logs/start_time.txt"
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
        "benchmarks/" + config["sub_dirs"]["trim_dir"] + "/trimmomatic_PE/{sample}.txt"
    shell:
        """
        trimmomatic PE \
            -threads {threads} \
            {input.reads} {output.R1_P} {output.R1_U} {output.R2_P} {output.R2_U} \
            ILLUMINACLIP:{params.adapters}:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:{params.qual} MINLEN:{params.minlen} \
            2> {log}
        """

rule phix_screen:
    message:
        """
        ** preprocessing **
        Removing phiX reads from {wildcards.sample}
        """
    input:
        R1 = config["sub_dirs"]["trim_dir"] + "/{sample}_1P.fastq.gz",
        R2 = config["sub_dirs"]["trim_dir"] + "/{sample}_2P.fastq.gz",
    output:
        R1 = config["sub_dirs"]["trim_dir"] + "/{sample}_1P.phiX.fastq.gz",
        R2 = config["sub_dirs"]["trim_dir"] + "/{sample}_2P.phiX.fastq.gz",
    log:
        "logs/phix_removal/{sample}.log"
    params:
        phix_genome = config["program_dir"] + config["phix_genome"]
    threads: 4
    shell:
        """
        bbduk.sh \
            in1={input.R1} \
            in2={input.R2} \
            outu1={output.R1} \
            outu2={output.R2} \
            threads={threads} \
            ref={params.phix_genome} \
            k=31 \
            hdist=1 \
            1>{log} 2>&1
        """


rule summarise_trimmomatic_log:
    input:
        expand("logs/trimmomatic_PE/{sample}.log", sample=config["samples"])
    output:
        "logs/trimmomatic_PE/trim_logs.summary"
    shell:
        """
        {config[program_dir]}/scripts/summarise_trimmomatic.py \
        -i {input} -o {output}
        """


rule porechop:
    message:
        """
        ** preprocessing **
        Trimming {wildcards.sample} MinION reads for adapters using Porechop
        """
    input:
        reads = getFastq
    output:
        config["sub_dirs"]["trim_dir"] + "/{sample}.porechop.fastq"
    params:
 
    threads: 8
    log:
        "logs/porechop/{sample}.log"
    benchmark:
        "benchmarks/porechop/{sample}.txt"
    shell:
        """
        porechop \
            -t {threads} \
            -i {input} \
            -o {output} \
            > {log}
        """
        
        
rule nanofilt:
    message:
        """
        ** preprocessing **
        Trimming {wildcards.sample} MinION reads for quality using NanoFilt
        """
    input:
        config["sub_dirs"]["trim_dir"] + "/{sample}.porechop.fastq",
    output:
        config["sub_dirs"]["trim_dir"] + "/{sample}.porechop.nanofilt.fastq",
    params:
        nanofilt_quality = config["nanofilt_quality"]
    threads: 1
    log:
        "logs/nanofilt/{sample}.log"
    benchmark:
        "benchmarks/nanofilt/{sample}.txt"
    shell:
        """
        NanoFilt \
            --quality {params.nanofilt_quality} \
            --logfile {log} \
            {input} \
            > {output}
        """

