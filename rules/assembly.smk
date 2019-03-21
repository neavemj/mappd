"""
Assembly rules

The inputs are:
    - trimmed reads in fastq format
The outputs are:
    - assembled contigs in fasta format
"""

rule spades:
    message:
        """
        Assembling trimmed reads with spades PE mode
        Using spades {config[data_type]} assembly
        This can be changed in the config file
        """
    input:
        R1_P = "01_trimmomatic/{sample}_1P.fastq.gz",
        R2_P = "01_trimmomatic/{sample}_2P.fastq.gz",
    output:
        directory("02_spades/{sample}")
    log:
        "logs/02_spades/{sample}.log"
    benchmark:
        "benchmarks/02_spades/{sample}.txt"
    threads: 8
    run:
        if config["data_type"] == "DNA":
            shell(
            """
            spades.py \
                -1 {input.R1_P} \
                -2 {input.R2_P} \
                -t {threads} \
                -m 8 \
                -o {output} > {log}
            """)
        elif config["data_type"] == "RNA":
            shell(
            """
            spades.py \
                -1 {input.R1_P} \
                -2 {input.R2_P} \
                -t {threads} \
                -m 8 \
                --rna \
                -o {output} > {log}
            """)
