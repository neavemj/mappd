"""
Host depletion rules

These rules will map reads to host genetic data
and split the reads into those that match and those that don't

Does a small assembly on a subset of the data
and blast the results to find host.
"""

rule subset_mRNA_reads:
    message:
        """
        ** host_depletion **
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

rule assemble_mRNA_subset:
    message:
        """
        ** host_depletion **
        Assembling small subset of {wildcards.sample} reads to identify host
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
        "benchmarks/" + config["sub_dirs"]["depletion_dir"] + "/spades_sub_assembly/{sample}.txt"
    threads: 16
    shell:
        """
        spades.py \
            -1 {input.R1} \
            -2 {input.R2} \
            -t {threads} \
            -o {params.out_dir} > {log}
        """

rule subset_subcontigs:
    message:
        """
        ** host_depletion **
        Gathering the 10 largest contigs from the sub-assembly
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_sub_assembly/contigs.fasta",
    output:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_largest_contigs.fasta",
    shell:
        # using the most abundant 20 contigs larger than 1000 bps for host identification
        # can experiment with these parameters if required
        """
        {config[program_dir]}/scripts/gather_contigs.py \
            -c {input} \
            -s 1000 \
            -n 20 \
            -o {output}
        """

rule blast_subcontigs:
    message:
        """
        ** host_depletion **
        Blasting the most abundant contigs from sub-assembly
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_largest_contigs.fasta",
    output:
        config["sub_dirs"]["depletion_dir"] + "/host/{sample}_largest_contigs.blastn",
    params:
        blast_nt = config["blast_nt"]
    threads: 16
    log:
        "logs/blast_sub_assembly/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["depletion_dir"] + "/blast_sub_assembly/{sample}.txt"
    shell:
        # in the output fmt, cols 6 and 7 need to be bitscore and taxid
        # for the scripts subset_blast.py and tally_organism_hits.py
        """
        blastn \
            -query {input} \
            -out {output} \
            -db {params.blast_nt} \
            -evalue 0.001 \
            -num_threads {threads} \
            -outfmt '6 \
                qseqid \
                sseqid \
                pident \
                length \
                evalue \
                bitscore \
                staxid \
                stitle'
        """

rule subset_subblast:
    message:
        """
        ** host_depletion **
        Retieving the 'best' hits for each {wildcards.sample} contig
        using the maximum bitscore
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

rule tally_abundant_subspecies:
    message:
        """
        ** host_depletion **
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
        {config[program_dir]}/scripts/tally_organism_hits.py \
            -b {input} \
            -t {output.wide} \
            -l {output.long}
        """

# TODO: If the host is identified as a virus or bacteria (maybe because the sample is a cell-culture)
# TODO: should stop the host-depletion and go direct to assembly.smk

rule associate_hostTaxid_genbank:
    message:
        """
        ** host_depletion **
        Retreiving genbank ids for the most abundant species
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/largest_contigs.blastn.tax.wide"
    output:
        config["sub_dirs"]["depletion_dir"] + "/host/host_nucl_nt.ids"
    params:
        # nt_to_taxids doesn't get as much info as including the wgs stuff
        # nt_to_taxids = config["nt_to_taxids"],
        gb_wgs_taxids = config["gb_wgs_taxids"],
        hosts_to_download = config["hosts_to_download"]
    benchmark:
        "benchmarks/" + config["sub_dirs"]["depletion_dir"] + "/grep_nucl_gb_ids/generic.txt"
    shell:
        # -w means only match whole words - don't want part of the taxid matching another one
        # cuts the first column (taxids), removes header,
        # and retrieves as many host taxids as required,
        """
        grep \
            -w \
            -f <(cut -f 1 {input} | tail -n +2 | head -n {params.hosts_to_download}) \
            {params.gb_wgs_taxids} \
            > {output}
        """

rule extract_host_nucl:
    message:
        """
        ** host_depletion **
        Extracting host nucleotide sequence from the nt database
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/host/host_nucl_nt.ids"
    output:
        config["sub_dirs"]["depletion_dir"] + "/host/host_nucl_nt.fasta"
    params:
        blast_nt = config["blast_nt"]
    log:
        "logs/extract_nucl_gb_fasta/accessions_not_found.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["depletion_dir"] + "/extract_nucl_gb_fasta/generic.txt"
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

rule bbmap_to_host:
    message:
        """
        ** host_depletion **
        Mapping mRNA {wildcards.sample} reads to host database
        """
    input:
        R1 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_mRNA_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_mRNA_2P.fastq",
        db = config["sub_dirs"]["depletion_dir"] + "/host/host_nucl_nt.fasta",
    output:
        R1 = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host_depleted_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host_depleted_2P.fastq",
    params:
        max_memory = config['bbmap_max_memory'],
    log:
        "logs/bbmap_host/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["depletion_dir"] + "/bbmap_host/{sample}.txt"
    threads: 16
    shell:
        """
        bbduk.sh \
            in1={input.R1} \
            in2={input.R2} \
            outu1={output.R1} \
            outu2={output.R2} \
            threads={threads} \
            {params.max_memory} \
            ref={input.db} \
            1>{log} 2>&1
        """


rule summarise_host_mapping:
    message:
        """
        ** host_depletion **
        Calculating number of reads mapped to host sequence
        """
    input:
        wide = config["sub_dirs"]["depletion_dir"] + "/host/largest_contigs.blastn.tax.wide",
        SSU_depleted = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_mRNA_1P.fastq",
        host_depleted = config["sub_dirs"]["depletion_dir"] + "/host/{sample}_host_depleted_1P.fastq",
    output:
        "logs/host_summary/{sample}_host_mapping_summary.tsv"
        #return(config["sub_dirs"]["depletion_dir"] + "/host/{}_host.blastn.abundance".format(wildcards.sample))
    shell:
        # first calculates how many reads in each file (divide by 2 gives all reads; divide by 4 to get pairs)
        # then uses these numbers to calculated actual LSU and SSU values
        # then gets taxid for the host from the wide format table above
        # then write this info to file
        """
        SSU_depleted_reads=$(expr $(cat {input.SSU_depleted} | wc -l) / 2)
        host_depleted_reads=$(expr $(cat {input.host_depleted} | wc -l) / 2)

        host_reads=$(expr $SSU_depleted_reads - $host_depleted_reads)

        taxid=$(cut -f 1 {input.wide} | tail -n +2 | head -n 1)

        echo -e \
            {wildcards.sample}'\t'${{taxid}}'\t'${{host_reads}} \
                > {output}
        """
