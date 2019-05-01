"""
Host depletion rules

These rules will map reads to host genetic data
and split the reads into those that match and those that don't
Using results from the rRNA mapping wasn't working particularly
well. Mainly because the SSU gene is too similar across species.
E.g. sheep and human are very similar, so the host species
gets mixed up.

Decided instead to do a small assembly on a subset of the data
and blast the results to find host.
"""

rule subset_mRNA_reads:
    message:
        """
        Taking a subset of 100,000 mRNA reads from {wildcards.sample}
        to be used for assembly and host identification
        """
    input:
        R1 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_mRNA_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_mRNA_2P.fastq",
    output:
        R1 = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_100k_mRNA_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_100k_mRNA_2P.fastq",
    shell:
        """
        head \
            -n 400000 \
            {input.R1} \
            > {output.R1} && \
        head \
            -n 400000 \
            {input.R2} \
            > {output.R2}
        """

rule assembling_mRNA_subset:
    message:
        """
        Assembling small subset {wildcards.sample} of reads to identify host
        """
    input:
        R1 = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_100k_mRNA_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_100k_mRNA_2P.fastq"
    output:
        out_fasta = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_sub_assembly/contigs.fasta",
    params:
        out_dir = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_sub_assembly"
    log:
        "logs/spades_sub_assembly/{sample}.log"
    benchmark:
        "benchmarks/spades_sub_assembly/{sample}.txt"
    threads: 16
    shell:
        """
        spades.py \
            -1 {input.R1} \
            -2 {input.R2} \
            -t {threads} \
            -o {params.out_dir} > {log}
        """

rule subset_contigs:
    message:
        """
        Gathering the 10 largest contigs from the sub-assembly
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_sub_assembly/contigs.fasta",
    output:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_largest_contigs.fasta",
    shell:
        # using the most abundant 10 contigs larger than 1000 bps for host identification
        # can experiment with this parameters if required
        """
        {config[program_dir]}/scripts/gather_contigs.py \
            -c {input} \
            -s 1000 \
            -n 10 \
            -o {output}
        """

rule blast_contigs:
    message:
        """
        Blasting the most abundant contigs from sub-assembly
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_largest_contigs.fasta",
    output:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_largest_contigs.blastn",
    params:
        blast_nt = config["blast_nt"]
    threads: 16
    shell:
        """
        blastn \
            -query {input} \
            -out {output} \
            -db {params.blast_nt} \
            -evalue 0.001 \
            -num_threads {threads} \
            -outfmt '6 qseqid sseqid pident length evalue bitscore staxid ssciname scomname salltitles'
        """

rule subset_blast:
    message:
        """
        Retieving the 'best' hits for each {wildcards.sample} contig
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_largest_contigs.blastn",
    output:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_largest_contigs.blastn.best_hits",
    shell:
        """
        {config[program_dir]}/scripts/subset_blast.py \
            -b {input} \
            -o {output}
        """

rule tally_abundant_species:
    message:
        """
        Calculating the most abundant species in the blast results for all samples
        """
    input:
        expand(config["sub_dirs"]["depletion_dir"] + "/host/{sample}_largest_contigs.blastn.best_hits", sample=config["samples"])
    output:
        # producing both wide and long format tables here
        # the wide will be used for the report, and the long for plotting in ggplot
        wide = config["sub_dirs"]["depletion_dir"] + "/host/largest_contigs.blastn.tax.wide",
        long = config["sub_dirs"]["depletion_dir"] + "/host/largest_contigs.blastn.tax.long"
    shell:
        """
        {config[program_dir]}/scripts/tally_abundant_hosts.py \
            -b {input} \
            -t {output.wide} \
            -l {output.long}
        """

rule plot_abundant_species:
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

rule associate_hostTaxid_genbank:
    message:
        """
        Retreiving genbank ids for the most abundant species
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/largest_contigs.blastn.tax.wide"
    output:
        config["sub_dirs"]["depletion_dir"] + "/host/host_nucl_gb.ids"
    params:
        acc_to_taxids = config["acc_to_taxids"],
        # need to do this to ensure only the string matches
        sed_pat = r"s/\(.*\)/\t\1\t/g",
        hosts_to_download = config["hosts_to_download"]
    benchmark:
        "benchmarks/grep_nucl_gb_ids/grep_nucl_gb_ids.txt"
    shell:
        # cuts the first column (taxids), removes header,
        # retrieves as many host taxids as required,
        # adds a tab before and after each taxid to ensure a
        # clean match
        """
        grep \
            -f <(cut -f 1 {input} | tail -n +2 | head -n {params.hosts_to_download} | sed "{params.sed_pat}") \
            {params.acc_to_taxids} \
            > {output}
        """

rule extract_host_nucl:
    message:
        """
        Extracting host nucleotide sequence from the nt database
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/host_nucl_gb.ids"
    output:
        config["sub_dirs"]["depletion_dir"] + "/host/host_nucl_gb.fasta"
    params:
        blast_nt = config["blast_nt"]
    log:
        "logs/extract_nucl_gb_fasta/accessions_not_found.log"
    benchmark:
        "benchmarks/extract_nucl_gb_fasta/extract_nucl_gb_fasta.txt"
    shell:
        # if this command doesn't find an accesion number (often)
        # it prints an error and returns an exit code of 1
        # this stops snakemake working
        # I want it to keep running even if some accessions weren't found
        # so the '|| true' bit ensures it returns a successful exit code
        """
        blastdbcmd \
            -db {params.blast_nt} \
            -entry_batch <(cut -f 2 {input}) \
            > {output} 2> {log} || true
        """

rule build_host_bowtiedb:
    message:
        """
        Building a bowtie2 database from host nucleotide sequence
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/host_nucl_gb.fasta"
    output:
        # bowtie2-build needs a basename for the database
        # usually I just give it the same name as the input
        # and it appends several *bt2 files
        # will trick snakemake by using this as an output even though
        # I won't use it in the shell command
        config["sub_dirs"]["depletion_dir"] + "/host/host_nucl_gb.fasta.1.bt2"
    shell:
        # use the same name for basename reference database
        """
        bowtie2-build \
            {input} \
            {input} > /dev/null
        """

rule bowtie_to_host:
    message:
        """
        Mapping {wildcards.sample} mRNA reads to host database
        """
    input:
        R1 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_mRNA_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_mRNA_2P.fastq",
        db_trick = config["sub_dirs"]["depletion_dir"] + "/host/host_nucl_gb.fasta.1.bt2"
    output:
        sam_fl = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host.sam"
    params:
        host_db = config["sub_dirs"]["depletion_dir"] + "/host/host_nucl_gb.fasta"
    log:
        "logs/bowtie_host/{sample}.log"
    benchmark:
        "benchmarks/bowtie_host/{sample}.txt"
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
        Converting {wildcards.sample} host sam file to bam
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host.sam"
    output:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host.bam"
    shell:
        """
        samtools view \
            -S -b \
            {input} > {output}
        """

rule host_get_unmapped:
    message:
        """
        Collecting {wildcards.sample} reads that did not map to host sequences
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host.bam"
    output:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host_depleted.bam"
    shell:
        # -f 13 should get reads where neither pair mapped (UNMAP & MUNMAP)
        """
        samtools view \
            -b \
            -f 13 \
            {input} > {output}
        """

rule host_sam_to_fastq:
    message:
        """
        Converting {wildcards.sample} host depleted sam file to fastq files
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host_depleted.bam"
    output:
        R1 = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host_depleted_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host_depleted_2P.fastq"
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

rule summarise_rRNA_host_mapping:
    message:
        """
        Summarising number of reads mapped to rRNA and host databases
        """
    input:
        lsu = expand("logs/bowtie_LSU/{sample}.log", sample=config["samples"]),
        ssu = expand("logs/bowtie_SSU/{sample}.log", sample=config["samples"]),
        host = expand("logs/bowtie_host/{sample}.log", sample=config["samples"])
    output:
        "logs/mapping_summary.tsv"
    shell:
        """
        {config[program_dir]}/scripts/summarise_mapping_results.py \
            -l {input.lsu} \
            -s {input.ssu} \
            -t {input.host} \
            -o {output}
        """

rule plot_mapping:
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


