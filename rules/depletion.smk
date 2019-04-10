"""
Depletion rules

These rules will map reads to the SILVA rRNA database
and split the reads into those that match and those that don't

"""

rule bowtie_to_rRNA:
    message:
        """
        Mapping cleaned reads to the SILVA rRNA database
        """
    input:
        R1_P = config["sub_dirs"]["trim_dir"] + "/{sample}_1P.fastq.gz",
        R2_P = config["sub_dirs"]["trim_dir"] + "/{sample}_2P.fastq.gz"
    output:
        sam_fl = config["sub_dirs"]["depletion_dir"] + "/{sample}.rRNA.sam"
    params:
        silva_db = config['silva_db'],
        # in bowtie the % symbol will be replaced with 1 or 2 depending on the pair
        mRNA_pairs = config["sub_dirs"]["depletion_dir"] + "/{sample}_mRNA_%P.fastq.gz"
    log:
        "logs/bowtie_rRNA/{sample}.log"
    benchmark:
        "benchmarks/bowtie_rRNA/{sample}.txt"
    threads:
        16
    shell:
        """
        bowtie2 \
            -x {params.silva_db} \
            -1 {input.R1_P} \
            -2 {input.R2_P} \
            -p {threads} \
            # write paired-end reads that fail to align concordantly (presumably not rRNA)
            --un-conc-gz {params.mRNA_pairs} \
            # suppress SAM records for unaligned reads to save space?
            --no-unal \
            -S {output.sam_fl}
        """
