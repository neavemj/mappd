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
        sam_fl = config["sub_dirs"]["depletion_dir"] + "/{sample}.LSU.sam"
    params:
        silva_LSU_db = config['silva_LSU_db'],
        # in bowtie the % symbol will be replaced with 1 or 2 depending on the pair
        # have to put the output files here because snakemake doesn't recognise the % symbol
        # this will be expanded in bowtie
        LSU_depleted_pairs = config["sub_dirs"]["depletion_dir"] + "/{sample}_LSU_depleted_%P.fastq.gz"
    log:
        "logs/bowtie_LSU/{sample}.log"
    benchmark:
        "benchmarks/bowtie_LSU/{sample}.txt"
    threads:
        16
    shell:
        # --un-conc-gz: write paired-end reads that fail to align concordantly (presumably not rRNA)
        # --no-unal: suppress SAM records for unaligned reads to save space
        # NOTE: need to use 2> because bowtie outputs to stderr
        """
        bowtie2 \
            -x {params.silva_LSU_db} \
            -1 {input.R1_P} \
            -2 {input.R2_P} \
            -p {threads} \
            --un-conc-gz {params.LSU_depleted_pairs} \
            --no-unal \
            -S {output.sam_fl} 2> {log}
        """

rule bowtie_to_SSU:
    message:
        """
        Mapping cleaned reads to the SILVA SSU rRNA database
        """
    input:
        # not really using this sam file but need to trick snakmake into running the LSU rule
        # can't use the actual fastq files as input because of the % symbol discussed above
        config["sub_dirs"]["depletion_dir"] + "/{sample}.LSU.sam",
    output:
        sam_fl = config["sub_dirs"]["depletion_dir"] + "/{sample}.SSU.sam",
    params:
        silva_SSU_db = config['silva_SSU_db'],
        # in bowtie the % symbol will be replaced with 1 or 2 depending on the pair
        # have to put the output files here because snakemake doesn't recognise the % symbol
        # this will be expanded in bowtie
        LSU_R1 = config["sub_dirs"]["depletion_dir"] + "/{sample}_LSU_depleted_1P.fastq.gz",
        LSU_R2 = config["sub_dirs"]["depletion_dir"] + "/{sample}_LSU_depleted_2P.fastq.gz",
        mRNA_pairs = config["sub_dirs"]["depletion_dir"] + "/{sample}_mRNA_%P.fastq.gz"
    log:
        "logs/bowtie_SSU/{sample}.log"
    benchmark:
        "benchmarks/bowtie_SSU/{sample}.txt"
    threads:
        16
    shell:
        # --un-conc-gz: write paired-end reads that fail to align concordantly (presumably not rRNA)
        # --no-unal: suppress SAM records for unaligned reads to save space
        """
        bowtie2 \
            -x {params.silva_SSU_db} \
            -1 {params.LSU_R1} \
            -2 {params.LSU_R2} \
            -p {threads} \
            --un-conc-gz {params.mRNA_pairs} \
            --no-unal \
            -S {output.sam_fl} 2> {log}
        """
