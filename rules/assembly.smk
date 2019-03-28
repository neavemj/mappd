"""
Assembly rules

The inputs are:
    - trimmed reads in fastq format
The outputs are:
    - assembled contigs in fasta format

Current assemblers include:
    - SPAdes in both RNA and DNA mode
    - Trinity for RNA-Seq only
"""

rule spades:
    message:
        """
        Assembling trimmed reads with spades PE mode
        Using spades {config[data_type]} assembly
        This can be changed in the config file
        """
    input:
        R1_P = config["sub_dirs"]["trim_dir"] + "/{sample}_1P.fastq.gz",
        R2_P = config["sub_dirs"]["trim_dir"] + "/{sample}_2P.fastq.gz",
    output:
        out_dir = directory(config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly"),
        out_trans = config["sub_dirs"]["assembly_dir"] + "/spades/{sample}/transcripts.fasta",
        out_gfa = config["sub_dirs"]["assembly_dir"] + "/spades/{sample}/assembly_graph_with_scaffolds.gfa"
    # have to do a params because spades --rna puts the graph file in the k-mer directory
    # however, normal spades puts it in the top level directory
    # I'll point to the file here so that for an --rna run, I can mv the file up a directory
    params:
        graph_fl = config["sub_dirs"]["assembly_dir"] + "/spades/{sample}/K73/assembly_graph_with_scaffolds.gfa"
    log:
        "logs/spades/{sample}.log"
    benchmark:
        "benchmarks/spades/{sample}.txt"
    threads: 16
    run:
        if config["data_type"] == "DNA":
            shell(
            """
            spades.py \
                -1 {input.R1_P} \
                -2 {input.R2_P} \
                -t {threads} \
                -m 8 \
                -o {output.out_dir} > {log}
            """)
        elif config["data_type"] == "RNA":
            shell(
            """
            spades.py \
                -1 {input.R1_P} \
                -2 {input.R2_P} \
                -t {threads} \
                -k 73 \
                -m 8 \
                --rna \
                -o {output.out_dir} > {log} &&
                mv {params.graph_fl} ..
            """)


rule trinity:
    message:
        """
        Assembling RNA-Seq reads with Trinity
        """
    input:
        R1_P = config["sub_dirs"]["trim_dir"] + "/{sample}_1P.fastq.gz",
        R2_P = config["sub_dirs"]["trim_dir"] + "/{sample}_2P.fastq.gz",
    output:
        directory(config["sub_dirs"]["assembly_dir"] + "/trinity/{sample}")
    log:
       "logs/trinity/{sample}.log"
    benchmark:
        "benchmarks/trinity/{sample}.txt"
    params:
        max_memory = "16G"
    threads: 16
    shell:
        """
        Trinity
            --seqType fq \
            --CPU {threads} \
            --max_memory {params.max_memory} \
            --full_cleanup \
            --left {input.R1_P} \
            --right {input.R2_P} \
            --output {output}
        """

rule subset_spades_bandage:
    input:
        config["sub_dirs"]["assembly_dir"] + "/spades/{sample}/assembly_graph_with_scaffolds.gfa"
    output:
        config["sub_dirs"]["assembly_dir"] + "/spades/{sample}/assembly_graph_10x.gfa"
    shell:
        """
        Bandage reduce \
            {input} \
            {output} \
            --scope depthrange \
            --mindepth 10 \
            --maxdepth 1000000
        """

rule draw_spades_bandage:
    input:
        config["sub_dirs"]["assembly_dir"] + "/spades/{sample}/assembly_graph_10x.gfa"
    output:
        config["sub_dirs"]["assembly_dir"] + "/spades/{sample}/assembly_graph_10x.png"
    shell:
        """
        Bandage image \
            {input} \
            {output}
        """
