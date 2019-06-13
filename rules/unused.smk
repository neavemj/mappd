
# bowtie rules no longer needed for mapping to host

rule build_host_bowtiedb:
    message:
        """
        ** host_depletion **
        Building a bowtie2 database from host nucleotide sequence
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/host_nucl_nt.fasta"
    output:
        # bowtie2-build needs a basename for the database
        # usually I just give it the same name as the input
        # and it appends several *bt2 files
        # will trick snakemake by using this as an output even though
        # I won't use it in the shell command
        config["sub_dirs"]["depletion_dir"] + "/host/host_nucl_nt.fasta.1.bt2"
    threads: 8
    shell:
        # use the same name for basename reference database
        """
        bowtie2-build \
            --threads {threads} \
            {input} \
            {input} > /dev/null
        """

rule bowtie_to_host:
    message:
        """
        ** host_depletion **
        Mapping {wildcards.sample} mRNA reads to host database
        """
    input:
        R1 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_mRNA_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_mRNA_2P.fastq",
        db_trick = config["sub_dirs"]["depletion_dir"] + "/host/host_nucl_nt.fasta.1.bt2"
    output:
        sam_fl = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host.sam"
    params:
        host_db = config["sub_dirs"]["depletion_dir"] + "/host/host_nucl_nt.fasta"
    log:
        "logs/bowtie_host/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["depletion_dir"] + "/bowtie_host/{sample}.txt"
    threads:
        16
    shell:
        """
        bowtie2 \
            -x {params.host_db} \
            -1 {input.R1} \
            -2 {input.R2} \
            -p {threads} \
            -S {output.sam_fl} 2> {log}
        """

rule host_sam_to_bam:
    message:
        """
        ** host_depletion **
        Converting {wildcards.sample} host sam file to bam
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host.sam"
    output:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host.bam"
    threads: 8
    shell:
        """
        samtools view \
            -@ {threads} \
            -S -b \
            {input} > {output}
        """

rule host_get_unmapped:
    message:
        """
        ** host_depletion **
        Collecting {wildcards.sample} reads that did not map to host sequences
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host.bam"
    output:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host_depleted.bam"
    threads: 8
    shell:
        # -f 13 should get reads where neither pair mapped (UNMAP & MUNMAP)
        """
        samtools view \
            -@ {threads} \
            -b \
            -f 13 \
            {input} > {output}
        """

rule host_sam_to_fastq:
    message:
        """
        ** host_depletion **
        Converting {wildcards.sample} host depleted sam file to fastq files
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host_depleted.bam"
    output:
        R1 = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host_depleted_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host_depleted_2P.fastq"
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


# this rule was copied over from the assembly.smk rules
# need to do this so I can add the host stats to the overall tax stats
rule host_mapping_stats:
    message:
        """
        ** host_depletion **
        Tallying statistics on {wildcards.sample} reads mapped to the host
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host.bam"
    output:
        sorted_bam = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host.sorted.bam",
        stats = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host.sorted.idxstats",
        # not including depth at this stage
        #depth = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host.sorted.depth",
    threads: 16
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
            > {output.stats}
        """

# required so that I can add the host reads onto my
# tax barcharts in the report
# had to modify the 'tally_organism_abundance.py' script here
# because I'm not adding reads per contig - rather per host taxid
rule summarise_host_abundance:
    message:
        """
        ** host_depletion **
        summarise abundance stats for the host species
        """
    input:
        wide = config["sub_dirs"]["depletion_dir"] + "/host/largest_contigs.blastn.tax.wide",
        stats = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host.sorted.idxstats",
        #depth = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host.sorted.depth",
        mapping = "logs/rRNA_mapping_summary.tsv",
    output:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host.blastn.abundance",
    shell:
        """
        {config[program_dir]}/scripts/tally_host_abundance.py \
            -w {input.wide} \
            -i {input.stats} \
            -m {input.mapping} \
            -o {output} \
        """






# rRNA depletion rules using bowtie / samtools
# replaced all this with single bbmap rules!

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

rule LSU_bam_to_fastq:
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

rule summarise_rRNA_mapping:
    message:
        """
        ** rRNA_depletion **
        Summarising number of reads mapped to rRNA databases
        """
    input:
        lsu = expand("logs/bowtie_LSU/{sample}.log", sample=config["samples"]),
        ssu = expand("logs/bowtie_SSU/{sample}.log", sample=config["samples"]),
    output:
        "logs/rRNA_mapping_summary.tsv"
    shell:
        """
        {config[program_dir]}/scripts/summarise_rRNA_mapping.py \
            -l {input.lsu} \
            -s {input.ssu} \
            -o {output}
        """

# can't get RNA simulation using rnftools to work properly

import rnftools

rnftools.mishmash.sample("simple_example",reads_in_tuple=2)

fa="GCA_003024735.1_ASM302473v1_cds_from_genomic.fna"

#contig_headers = [contig.strip().lstrip(">").split(" ")[0] for contig in open(fa) if contig.startswith(">")]

reads_required = 10000

# this doesn't work for some reason - weird 'pickling' error

for i in [2, 3, 7]:
    rnftools.mishmash.ArtIllumina(
    	fasta=fa,
    	#number_of_read_tuples=10000,
        sequences=[i],
        coverage=10,
    	read_length_1=150,
    	read_length_2=150,
    )

include: rnftools.include()
rule: input: rnftools.input()


# rRNA depletion rules no longer required

rRNA_type = ["LSU", "SSU"]

rule associate_genbank_to_taxids:
    message:
        """
        ** rRNA_depletion **
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
        "benchmarks/" + config["sub_dirs"]["depletion_dir"] + "/grep_taxids_{rRNA_type}/{sample}.txt"
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
        ** rRNA_depletion **
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
        "benchmarks/" + config["sub_dirs"]["depletion_dir"] + "/grep_silva{rRNA_type}_taxstrings/{sample}.txt"
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
        ** rRNA_depletion **
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


# unused host depletion rules

rule plot_abundant_subspecies:
    input:
        long = config["sub_dirs"]["depletion_dir"] + "/host/largest_contigs.blastn.tax.long"
    output:
        pdf = config["sub_dirs"]["depletion_dir"] + "/host/largest_contigs.blastn.tax.pdf",
        png = config["sub_dirs"]["depletion_dir"] + "/host/largest_contigs.blastn.tax.png"
    shell:
        """
        Rscript {config[program_dir]}/scripts/plot_host_stats.R \
        {input} {output.pdf} {output.png}
        """

rule plot_host_mapping:
    input:
        "logs/mapping_summary.tsv"
    output:
        pdf = "logs/mapping_summary.pdf",
        png = "logs/mapping_summary.png"
    shell:
        """
        Rscript {config[program_dir]}/scripts/plot_mapping.R \
            {input} {output.pdf} {output.png}
        """

# unused assembly rules

rule subset_spades_bandage:
    input:
        config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/assembly_graph_with_scaffolds.gfa"
    output:
        config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/assembly_graph_10x.gfa"
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
        config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/assembly_graph_10x.gfa"
    output:
        config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/assembly_graph_10x.png"
    shell:
        """
        Bandage image \
            {input} \
            {output}
        """

rule draw_trinity_bandage:
    input:
        config["sub_dirs"]["assembly_dir"] + "/trinity/{sample}_trinity/{sample}_trinity.Trinity.fasta"
    output:
        config["sub_dirs"]["assembly_dir"] + "/trinity/{sample}_trinity/{sample}_trinity.Trinity.png"
    shell:
        """
        Bandage image \
            {input} \
            {output}
        """
