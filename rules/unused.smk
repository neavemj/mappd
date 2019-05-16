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
