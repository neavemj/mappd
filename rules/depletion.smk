"""
Depletion rules

These rules will map reads to the SILVA LSU and SSU rRNA databases
and split the reads into those that match and those that don't

"""

rule bowtie_to_LSU:
    message:
        """
        Mapping cleaned reads to the SILVA LSU rRNA database
        """
    input:
        R1_P = config["sub_dirs"]["trim_dir"] + "/{sample}_1P.fastq.gz",
        R2_P = config["sub_dirs"]["trim_dir"] + "/{sample}_2P.fastq.gz"
    output:
        sam_fl = config["sub_dirs"]["depletion_dir"] + "/{sample}_LSU.sam"
    params:
        silva_LSU_db = config['silva_LSU_db']
    log:
        "logs/bowtie_LSU/{sample}.log"
    benchmark:
        "benchmarks/bowtie_LSU/{sample}.txt"
    threads:
        16
    shell:
        # NOTE: need to use 2> because bowtie outputs to stderr
        """
        bowtie2 \
            -x {params.silva_LSU_db} \
            -1 {input.R1_P} \
            -2 {input.R2_P} \
            -p {threads} \
            -S {output.sam_fl} 2> {log}
        """

rule LSU_get_unmapped:
    message:
        """
        Collecting reads that did not map to the LSU database
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/{sample}_LSU.sam"
    output:
        config["sub_dirs"]["depletion_dir"] + "/{sample}_LSU_depleted.sam"
    shell:
        # -f 13 should get reads where neither pair mapped (UNMAP & MUNMAP)
        """
        samtools view \
            -f 13 \
            {input} > {output}
        """

rule LSU_sam_to_fastq:
    message:
        """
        Converting LSU depleted sam file to fastq files
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/{sample}_LSU_depleted.sam"
    output:
        R1 = config["sub_dirs"]["depletion_dir"] + "/{sample}_LSU_depleted_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/{sample}_LSU_depleted_2P.fastq"
    shell:
    # the dev null bits put discard unpaired reads
    # the -F bit ensures the mates are paired
        """
        samtools fastq \
            -1 {output.R1} \
            -2 {output.R2} \
            -0 /dev/null \
            -s /dev/null \
            -n \
            -F 0x900 \
            {input} > /dev/null
        """

rule bowtie_to_SSU:
    message:
        """
        Mapping LSU-depleted reads to the SILVA SSU rRNA database
        """
    input:
        R1 = config["sub_dirs"]["depletion_dir"] + "/{sample}_LSU_depleted_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/{sample}_LSU_depleted_2P.fastq"
    output:
        sam_fl = config["sub_dirs"]["depletion_dir"] + "/{sample}_SSU.sam"
    params:
        silva_SSU_db = config['silva_SSU_db'],
    log:
        "logs/bowtie_SSU/{sample}.log"
    benchmark:
        "benchmarks/bowtie_SSU/{sample}.txt"
    threads:
        16
    shell:
        """
        bowtie2 \
            -x {params.silva_SSU_db} \
            -1 {input.R1} \
            -2 {input.R2} \
            -p {threads} \
            -S {output.sam_fl} 2> {log}
        """

rule SSU_get_unmapped:
    message:
        """
        Collecting reads that did not map to either the LSU or SSU database
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/{sample}_SSU.sam"
    output:
        config["sub_dirs"]["depletion_dir"] + "/{sample}_rRNA_depleted.sam"
    shell:
        # -f 13 should get reads where neither pair mapped (UNMAP & MUNMAP)
        """
        samtools view \
            -f 13 \
            {input} > {output}
        """

rule mRNA_sam_to_fastq:
    message:
        """
        Converting rRNA depleted sam file to mRNA fastq files
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/{sample}_rRNA_depleted.sam"
    output:
        R1 = config["sub_dirs"]["depletion_dir"] + "/{sample}_mRNA_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/{sample}_mRNA_2P.fastq"
    shell:
    # the dev null bits put discard unpaired reads
    # the -F bit ensures the mates are paired
        """
        samtools fastq \
            -1 {output.R1} \
            -2 {output.R2} \
            -0 /dev/null \
            -s /dev/null \
            -n \
            -F 0x900 \
            {input} > /dev/null
        """
