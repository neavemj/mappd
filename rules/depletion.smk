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
        sam_fl = config["sub_dirs"]["depletion_dir"] + "/{sample}.LSU.sam",
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
        """
        bowtie2 \
            -x {params.silva_LSU_db} \
            -1 {input.R1_P} \
            -2 {input.R2_P} \
            -p {threads} \
            --un-conc-gz {params.LSU_depleted_pairs} \
            --no-unal \
            -S {output.sam_fl} > {log}
        """

