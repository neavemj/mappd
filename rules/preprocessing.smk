"""
These rules will trim sequencing reads

The inputs are:
    - raw Illumina sequencing reads
The outputs are:
    - cleaned Illumina sequencing reads
"""

import os, sys
from os.path import join
import json
import glob

config: "config.yaml"

# record start time of the run

rule record_start:
    output:
        "logs/start_time.txt"
    shell:
        """
        echo -e "start time\t"$(date) > {output}
        """

# 2025-07-07 updated to get sample names from directly from raw file directory
RAW_DIR = config["raw_dir"]

# need to differentiate between Nextseq or MiSeq data
# due to slight differences in the file naming
# using a regex so that the run sample number is not captured "i.e, S2"
if config["illumina_machine"] == "nextseq":
    SAMPLES, = glob_wildcards(join(RAW_DIR, "{sample,[^_]+(?:-[^_]+)*}_S[0-9]+_R1_001.fastq.gz"))
elif config["illumina_machine"] == "miseq":
    SAMPLES, = glob_wildcards(join(RAW_DIR, "{sample,[^_]+(?:-[^_]+)*}_S[0-9]+_L001_R1_001.fastq.gz"))


# need a function here to differentiate nextseq from miseq dataset
# also need to add back in the S* that I omitted earlier
def get_reads(wildcards=None, sample=None):
    # Safely resolve sample name
    if sample is not None:
        sample_prefix = sample
    elif hasattr(wildcards, "sample"):
        sample_prefix = wildcards.sample
    elif isinstance(wildcards, str):
        sample_prefix = wildcards
    else:
        raise ValueError("get_reads() must be called with a 'sample' or 'wildcards.sample'")

    # Match files based on machine type
    if config["illumina_machine"] == "nextseq":
        r1 = glob.glob(join(RAW_DIR, f"{sample_prefix}_S*_R1_001.fastq.gz"))[0]
        r2 = glob.glob(join(RAW_DIR, f"{sample_prefix}_S*_R2_001.fastq.gz"))[0]
    elif config["illumina_machine"] == "miseq":
        r1 = glob.glob(join(RAW_DIR, f"{sample_prefix}_S*_L001_R1_001.fastq.gz"))[0]
        r2 = glob.glob(join(RAW_DIR, f"{sample_prefix}_S*_L001_R2_001.fastq.gz"))[0]
    else:
        raise ValueError(f"Unknown illumina_machine: {config['illumina_machine']}")
        
    return [r1, r2]
    
# want to create a sample config file for the reporting of samples / paths
rule create_sample_config:
    output:
        "logs/sample_config.json"
    run:
        sample_dict = {
            sample: get_reads(sample)
            for sample in SAMPLES
        }
        with open(output[0], "w") as f:
            json.dump(sample_dict, f, indent=2)

# TODO: might need a trimmomatic SE mode
rule trimmomatic_PE:
    message:
        """
        ** preprocessing **
        Trimming {wildcards.sample} for quality and Illumina adapters using Trimmomatic
        """
    input:
        reads = get_reads,
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
        lambda wildcards: expand("logs/trimmomatic_PE/{sample}.log", sample=SAMPLES)
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
        reads = get_reads
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

