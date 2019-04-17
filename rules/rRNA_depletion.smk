"""
Ribosomal RNA depletion rules

These rules will map reads to the SILVA LSU and SSU rRNA databases
and split the reads into those that match and those that don't
It will also guess the host species based on rRNA matches

"""

rule bowtie_to_LSU:
    message:
        """
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

rule LSU_sam_to_bam:
    message:
        """
        Converting {wildcards.sample} LSU sam file to bam
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU.sam"
    output:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU.bam"
    shell:
        """
        samtools view \
            -S -b \
            {input} > {output}
        """

rule LSU_stats:
    message:
        """
        Tallying statistics on {wildcards.sample} reads mapped to the LSU database
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU.bam"
    output:
        sorted_bam = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU.sorted.bam",
        stats = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU.idxstats"
    shell:
        # this will sort > index > idxstats >
        # sort by number of mapped reads > only output contigs with at least 3 read mapped
        # i.e. more than a singleton pair
        """
        samtools sort \
            {input} > {output.sorted_bam} && \
        samtools index \
            {output.sorted_bam} && \
        samtools idxstats \
            {output.sorted_bam} | \
            sort -nrk 3 | \
            awk '$3 > 2' > {output.stats}
        """

rule LSU_get_unmapped:
    message:
        """
        Collecting {wildcards.sample} reads that did not map to the LSU database
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU.bam"
    output:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU_depleted.bam"
    shell:
        # -f 13 should get reads where neither pair mapped (UNMAP & MUNMAP)
        # will turn into a bam file to save space
        """
        samtools view \
            -f 13 \
            {input} > {output}
        """

rule LSU_sam_to_fastq:
    message:
        """
        Converting {wildcards.sample} LSU depleted sam file to fastq files
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU_depleted.bam"
    output:
        R1 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU_depleted_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU_depleted_2P.fastq"
    shell:
    # the dev null bit discards unpaired reads
    # the -F bit ensures the mates are paired
        """
        samtools fastq \
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

rule SSU_sam_to_bam:
    message:
        """
        Converting {wildcards.sample} SSU sam file to bam
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_SSU.sam"
    output:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_SSU.bam"
    shell:
        """
        samtools view \
            -S -b \
            {input} > {output}
        """

rule SSU_stats:
    message:
        """
        Tallying {wildcards.sample} statistics on reads mapped to the SSU database
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_SSU.bam"
    output:
        sorted_bam = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_SSU.sorted.bam",
        stats = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_SSU.idxstats"
    shell:
        # this will sort > index > idxstats >
        # sort by number of mapped reads > only output contigs with at least 3 read mapped
        # i.e. more than a singleton pair
        """
        samtools sort \
            {input} > {output.sorted_bam} && \
        samtools index \
            {output.sorted_bam} && \
        samtools idxstats \
            {output.sorted_bam} | \
            sort -nrk 3 | \
            awk '$3 > 2' > {output.stats}
        """

rule SSU_get_unmapped:
    message:
        """
        Collecting {wildcards.sample} reads that did not map to either the LSU or SSU database
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_SSU.bam"
    output:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_rRNA_depleted.bam"
    shell:
        # -f 13 should get reads where neither pair mapped (UNMAP & MUNMAP)
        """
        samtools view \
            -f 13 \
            -b \
            {input} > {output}
        """

rule mRNA_sam_to_fastq:
    message:
        """
        Converting {wildcards.sample} rRNA depleted sam file to mRNA fastq files
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_rRNA_depleted.bam"
    output:
        R1 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_mRNA_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_mRNA_2P.fastq"
    shell:
    # the dev null bit discards unpaired reads
    # the -F bit ensures the mates are paired
        """
        samtools fastq \
            -1 {output.R1} \
            -2 {output.R2} \
            -0 /dev/null \
            -s /dev/null \
            -n \
            -F 0x900 \
            {input} 2> /dev/null
        """

rRNA_type = ["LSU", "SSU"]

rule associate_genbank_to_taxids:
    message:
        """
        Associating {wildcards.sample} genbank accessions from the SILVA database to NCBI taxids
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_{rRNA_type}.idxstats"
    output:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_{rRNA_type}.idxstats.taxid"
    params:
        # for some reason, I have to define the sed pattern here, then pass it in
        # otherwise it is weirdly expanded in the actual shell command
        sed_pat = r"s/\(.*\)/^\1\t/g",
        acc_to_taxids = config["acc_to_taxids"]
    benchmark:
        "benchmarks/grep_taxids/{sample}_{rRNA_type}.txt"
    shell:
        # takes the idxstats output, cuts the first column (genbank ID)
        # ensures only uniq entries, then adds '^' to the start of the pattern,
        # and \t to the end. This ensures a complete match.
        # then greps for this pattern in the acc_to_taxids file
        """
        grep \
            -f <(cut -f 1 -d "." {input} | sort | uniq | sed "{params.sed_pat}") \
            {params.acc_to_taxids} \
            > {output}
        """

rule associate_genbank_to_silvaids:
    message:
        """
        Associating {wildcards.sample} genbank accessions from the SILVA database to their taxonomy string
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_{rRNA_type}.idxstats"
    output:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_{rRNA_type}.idxstats.taxstring"
    params:
        # for some reason, I have to define the sed pattern here, then pass it in
        # otherwise it is weirdly expanded in the actual shell command
        # pat puts > back and stops at first ".". This matches silva format
        sed_pat = r"s/\(.*\)/>\1\./g",
        # created these files by simply grepping out the header line from the fasta files
        silva_LSU_db = config['silva_LSU_db'] + ".taxstring",
        silva_SSU_db = config['silva_SSU_db'] + ".taxstring",
    benchmark:
        "benchmarks/grep_silva_taxstrings/{sample}_{rRNA_type}.txt"
    run:
        # doing an if else statement here to use the correct database
        if wildcards.rRNA_type == "LSU":
            # takes the idxstats output, cuts the first column (genbank ID)
            # ensures only uniq entries, then adds '^' to the start of the pattern,
            # and \t to the end. This ensures a complete match.
            # then greps for this pattern in the acc_to_taxids file
            # -f <(cut -f 1 -d "." {input} | sort | uniq | sed "s/\(.*\)/^\1\t/g") \
            shell(
            """
            grep \
                -f <(cut -f 1 -d "." {input} | sort | uniq | sed "{params.sed_pat}") \
                {params.silva_LSU_db} \
                > {output}
            """)
        elif wildcards.rRNA_type == "SSU":
            shell(
            """
            grep \
                -f <(cut -f 1 -d "." {input} | sort | uniq | sed "{params.sed_pat}") \
                {params.silva_SSU_db} \
                > {output}
            """)

rule add_taxonomy_to_idxstats:
    message:
        """
        Adding taxid and SILVA taxonomy string to the idxstats file
        """
    input:
        idxstats = expand(config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_{rRNA_type}.idxstats", sample=config["samples"],
            rRNA_type=rRNA_type),
        taxstring = expand(config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_{rRNA_type}.idxstats.taxstring", sample=config["samples"],
            rRNA_type=rRNA_type),
        taxid = expand(config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_{rRNA_type}.idxstats.taxid", sample=config["samples"],
            rRNA_type=rRNA_type)
    output:
        full_table = config["sub_dirs"]["depletion_dir"] + "/rRNA/idxstats.summary"
    log:
        "logs/add_taxidTaxstring/taxid_not_found.log"
    shell:
        """
        {config[program_dir]}/scripts/summarise_rRNAstats.py \
            -i {input.idxstats} \
            -t {input.taxstring} \
            -d {input.taxid} \
            -o {output.full_table} > {log}
        """

rule plot_idxstats:
    # decided to just plot and output tables for the SSU results in the report
    # most people sequence the 18S gene so this should be more accurate
    # the LSU reads will still be removed, just not plotted
    # this could be easily changed in the below R script
    input:
        full_table = config["sub_dirs"]["depletion_dir"] + "/rRNA/idxstats.summary",
    output:
        pdf = config["sub_dirs"]["depletion_dir"] + "/rRNA/SSU.idxstats.summary.pdf",
        png = config["sub_dirs"]["depletion_dir"] + "/rRNA/SSU.idxstats.summary.png",
        tsv = config["sub_dirs"]["depletion_dir"] + "/rRNA/SSU.idxstats.summary.tsv"
    shell:
        # couldn't figure out in R how to add a simple space after the colon!
        # just do it in sed
        """
        Rscript {config[program_dir]}/scripts/plot_rRNA_idxstats.R \
        {input} {output.pdf} {output.png} {output.tsv} && \
        sed -i "s/;/; /g" {output.tsv}
        """



