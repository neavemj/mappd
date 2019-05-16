"""
Ribosomal RNA depletion rules

These rules will map reads to the SILVA LSU and SSU rRNA databases
and split the reads into those that match and those that don't
It will also guess the host species based on rRNA matches

"""

rule bowtie_to_LSU:
    message:
        """
        ** rRNA_depletion **
        Mapping cleaned {wildcards.sample} reads to the SILVA LSU rRNA database
        """
    input:
        R1_P = config["sub_dirs"]["trim_dir"] + "/{sample}_1P.fastq.gz",
        R2_P = config["sub_dirs"]["trim_dir"] + "/{sample}_2P.fastq.gz"
    output:
        # can mark these large sam files with temp() and they will be
        # deleted when no other rules need them anymore
        sam_fl = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU.sam"
    params:
        silva_LSU_db = config['silva_LSU_db']
    log:
        "logs/bowtie_LSU/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["depletion_dir"] + "/bowtie_LSU/{sample}.txt"
    threads: 16
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

rule LSU_sam_to_bam:
    message:
        """
        ** rRNA_depletion **
        Converting {wildcards.sample} LSU sam file to bam
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU.sam"
    output:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU.bam"
    threads: 8
    shell:
        """
        samtools view \
            -@ {threads} \
            -S -b \
            {input} > {output}
        """

rule LSU_stats:
    message:
        """
        ** rRNA_depletion **
        Tallying statistics on {wildcards.sample} reads mapped to the LSU database
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU.bam"
    output:
        sorted_bam = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU.sorted.bam",
        stats = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU.idxstats"
    threads: 8
    shell:
        # this will sort > index > idxstats >
        # sort by number of mapped reads > only output contigs with at least 3 read mapped
        # i.e. more than a singleton pair
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
            sort -nrk 3 | \
            awk '$3 > 2' > {output.stats}
        """

rule LSU_get_unmapped:
    message:
        """
        ** rRNA_depletion **
        Collecting {wildcards.sample} reads that did not map to the LSU database
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU.bam"
    output:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU_depleted.bam"
    threads: 8
    shell:
        # -f 13 should get reads where neither pair mapped (UNMAP & MUNMAP)
        # will turn into a bam file to save space
        """
        samtools view \
            -@ {threads} \
            -f 13 \
            {input} > {output}
        """

rule LSU_sam_to_fastq:
    message:
        """
        ** rRNA_depletion **
        Converting {wildcards.sample} LSU depleted sam file to fastq files
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU_depleted.bam"
    output:
        R1 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU_depleted_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU_depleted_2P.fastq"
    threads: 8
    shell:
    # the dev null bit discards unpaired reads
    # the -F bit ensures the mates are paired
        """
        samtools fastq \
            -@ {threads} \
            -1 {output.R1} \
            -2 {output.R2} \
            -0 /dev/null \
            -s /dev/null \
            -n \
            -F 0x900 \
            {input} 2> /dev/null
        """

rule bowtie_to_SSU:
    message:
        """
        ** rRNA_depletion **
        Mapping {wildcards.sample} LSU-depleted reads to the SILVA SSU rRNA database
        """
    input:
        R1 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU_depleted_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU_depleted_2P.fastq"
    output:
        sam_fl = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_SSU.sam"
    params:
        silva_SSU_db = config['silva_SSU_db']
    log:
        "logs/bowtie_SSU/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["depletion_dir"] + "/bowtie_SSU/{sample}.txt"
    threads: 16
    shell:
        """
        bowtie2 \
            -x {params.silva_SSU_db} \
            -1 {input.R1} \
            -2 {input.R2} \
            -p {threads} \
            -S {output.sam_fl} 2> {log}
        """

rule SSU_sam_to_bam:
    message:
        """
        ** rRNA_depletion **
        Converting {wildcards.sample} SSU sam file to bam
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_SSU.sam"
    output:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_SSU.bam"
    threads: 8
    shell:
        """
        samtools view \
            -@ {threads} \
            -S -b \
            {input} > {output}
        """

rule SSU_stats:
    message:
        """
        ** rRNA_depletion **
        Tallying {wildcards.sample} statistics on reads mapped to the SSU database
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_SSU.bam"
    output:
        sorted_bam = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_SSU.sorted.bam",
        stats = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_SSU.idxstats"
    threads: 8
    shell:
        # this will sort > index > idxstats >
        # sort by number of mapped reads > only output contigs with at least 3 read mapped
        # i.e. more than a singleton pair
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
            sort -nrk 3 | \
            awk '$3 > 2' > {output.stats}
        """

rule SSU_get_unmapped:
    message:
        """
        ** rRNA_depletion **
        Collecting {wildcards.sample} reads that did not map to either the LSU or SSU database
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_SSU.bam"
    output:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_rRNA_depleted.bam"
    threads: 8
    shell:
        # -f 13 should get reads where neither pair mapped (UNMAP & MUNMAP)
        """
        samtools view \
            -@ {threads} \
            -f 13 \
            -b \
            {input} > {output}
        """

rule mRNA_sam_to_fastq:
    message:
        """
        ** rRNA_depletion **
        Converting {wildcards.sample} rRNA depleted sam file to mRNA fastq files
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_rRNA_depleted.bam"
    output:
        R1 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_mRNA_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_mRNA_2P.fastq"
    threads: 8
    shell:
    # the dev null bit discards unpaired reads
    # the -F bit ensures the mates are paired
        """
        samtools fastq \
            -@ {threads} \
            -1 {output.R1} \
            -2 {output.R2} \
            -0 /dev/null \
            -s /dev/null \
            -n \
            -F 0x900 \
            {input} 2> /dev/null
        """

