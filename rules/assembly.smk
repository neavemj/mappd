"""
Assembly rules

The inputs are:
    - trimmed reads in fastq format
The outputs are:
    - assembled contigs in fasta format

Current assemblers include:
    - SPAdes in RNA mode
    - Trinity for RNA-Seq only
"""

# first need to decide which reads to assemble based on the config file
# options include rRNA depletion or host depletion
# the host depletion happens last so the following elif order should work
def get_reads(wildcards):
    # if host_depletion is true then rRNA_depletion will also be true because
    # the host_depletion rules use rRNA depleted fastq files
    if config["host_depletion"]:
        return([
            config["sub_dirs"]["depletion_dir"] + "/host/{}_host_depleted_1P.fastq".format(wildcards.sample),
            config["sub_dirs"]["depletion_dir"] + "/host/{}_host_depleted_2P.fastq".format(wildcards.sample)
        ])
    # if host_depletion is not true, can still just do the rRNA_depletion bit
    elif config["rRNA_depletion"]:
        return([
            config["sub_dirs"]["depletion_dir"] + "/rRNA/{}_mRNA_1P.fastq".format(wildcards.sample),
            config["sub_dirs"]["depletion_dir"] + "/rRNA/{}_mRNA_2P.fastq".format(wildcards.sample)
        ])
    # if neither host_depletion or rRNA_depletion are true, just returned trimmed reads
    # NOTE: At this stage, rRNA_depletion will be run regardless due to downstream requirements
    # need to modify overall_abundance plot (required by report) to make this work
    else:
        return([
            config["sub_dirs"]["trim_dir"] + "/{}_1P.fastq.gz".format(wildcards.sample),
            config["sub_dirs"]["trim_dir"] + "/{}_2P.fastq.gz".format(wildcards.sample)
        ])


rule spades:
    message:
        """
        ** assembly **
        Assembling {wildcards.sample} mRNA reads with spades PE mode
        Using spades RNA mode
        """
    input:
        get_reads
    output:
        out_trans = config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts.fasta",
        out_gfa = config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/assembly_graph_with_scaffolds.gfa"
    # have to do a params because spades --rna puts the graph file in the k-mer directory
    # however, normal spades puts it in the top level directory
    # I'll point to the file here so that for an --rna run, I can mv the file up a directory
    params:
        graph_fl = config["sub_dirs"]["assembly_dir"] + "/spades/{"
            "sample}_assembly/K73/assembly_graph_with_scaffolds.gfa",
        out_dir = config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly",
        max_memory = 16
    log:
        "logs/spades/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["assembly_dir"] + "/spades/{sample}.txt"
    threads: 16
    shell:
        """
        spades.py \
            -1 {input[0]} \
            -2 {input[1]} \
            -t {threads} \
            -k 73 \
            -m {params.max_memory} \
            --rna \
            -o {params.out_dir} > {log} &&
            mv {params.graph_fl} {params.out_dir}
        """

rule trinity:
    message:
        """
        ** assembly **
        Assembling RNA-Seq reads with Trinity
        """
    input:
        R1 = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host_depleted_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host_depleted_2P.fastq"
    output:
        config["sub_dirs"]["assembly_dir"] + "/trinity/{sample}_trinity/{sample}_trinity.Trinity.fasta"
    log:
       "logs/trinity/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["depletion_dir"] + "/trinity/{sample}.txt"
    params:
        out_dir = config["sub_dirs"]["assembly_dir"] + "/trinity/{sample}_trinity/{sample}_trinity",
        max_memory = "16G"
    threads: 16
    shell:
        """
        Trinity \
            --seqType fq \
            --CPU {threads} \
            --max_memory {params.max_memory} \
            --full_cleanup \
            --left {input.R1} \
            --right {input.R2} \
            --output {params.out_dir} > {log}
        """

rule subset_spades_contigs:
    message:
        """
        ** assembly **
        Removing small contigs from the assembly
        """
    input:
        config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts.fasta"
    output:
        config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts_subset.fasta"
    params:
        min_contig_size = config["min_contig_size"]
    shell:
        # extract contigs larger than 500 bps for annotation
        # stop at 100,000 contigs? Presumably there won't be many more that this?
        """
        {config[program_dir]}/scripts/gather_contigs.py \
            -c {input} \
            -s {params.min_contig_size} \
            -n 100000 \
            -o {output}
        """

rule build_spades_bowtiedb:
    message:
        """
        ** assembly **
        Building a bowtie2 database for the SPAdes assembly
        """
    input:
        config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts_subset.fasta"
    output:
        # bowtie2-build needs a basename for the database
        # usually I just give it the same name as the input
        # and it appends several *bt2 files
        # will trick snakemake by using this as an output even though
        # I won't use it in the shell command
        config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts_subset.fasta.1.bt2"
    threads: 8
    shell:
        # use the same name for basename reference database
        """
        bowtie2-build \
            --threads {threads} \
            {input} \
            {input} > /dev/null
        """

rule bowtie_to_spades_assembly:
    message:
        """
        ** assembly **
        Mapping {wildcards.sample} reads back to SPAdes assembly
        to get abundance estimates
        """
    input:
        # map either host_depleted reads or raw reads back to assembly
        # depending on config file requirements
        reads = get_reads,
        db_trick = config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts_subset.fasta.1.bt2"
    output:
        sam_fl = config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts_subset.sam"
    params:
        assembly_db = config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts_subset.fasta"
    log:
        "logs/bowtie_spades_assembly/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["assembly_dir"] + "/bowtie_spades_assembly/{sample}.txt"
    threads: 16
    shell:
        """
        bowtie2 \
            -x {params.assembly_db} \
            -1 {input.reads[0]} \
            -2 {input.reads[1]} \
            -p {threads} \
            -S {output.sam_fl} 2> {log}
        """

rule spades_sam_to_bam:
    message:
        """
        ** assembly **
        Converting {wildcards.sample} SPAdes sam file to bam
        """
    input:
        config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts_subset.sam"
    output:
        config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts_subset.bam"
    threads: 8
    shell:
        """
        samtools view \
            -@ {threads} \
            -S -b \
            {input} > {output}
        """

rule spades_mapping_stats:
    message:
        """
        ** assembly **
        Tallying statistics on {wildcards.sample} reads mapped to the SPAdes assembly
        """
    input:
        config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts_subset.bam"
    output:
        sorted_bam = config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts_subset.sorted.bam",
        stats = config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts_subset.sorted.idxstats",
        depth = config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts_subset.sorted.depth",
    threads: 8
    shell:
        # this will sort > index > idxstats > sort by most mapped reads
        """
        samtools sort \
            -@ {threads} \
            {input} > {output.sorted_bam} && \
        samtools index \
            -@ {threads} \
            {output.sorted_bam} && \
        samtools idxstats \
            -@ {threads} \
            {output.sorted_bam} | \
            sort -nrk 3 \
            > {output.stats} && \
        samtools depth \
            {output.sorted_bam} \
            > {output.depth}
        """

