"""
Annotation rules

Takes contigs assembled by Spades or Trinity
and does blast / diamond searches of the contigs
againt NCBI databases.

Options:
- blastx / diamondx => nr database (10^-5)
- blastx => viral protein RefSeq (10^-5)
    - then compare hits to complete nr database to
      ensure the match is still the best
"""


rule diamond_nr:
    message:
        """
        ** annotation **
        Using Diamond blastx to compare {wildcards.sample} contigs to the nr database
        using an evalue of {params.diamond_nr_evalue}
        """
    input:
        config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts_subset.fasta"
    output:
        config["sub_dirs"]["annotation_dir"] + "/diamond/{sample}.diamond_blastx"
    params:
        diamond_nr_db = config["diamond_nr"],
        diamond_nr_evalue = config["diamond_nr_evalue"]
    log:
       "logs/diamond/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["annotation_dir"] + "/diamond/{sample}.txt"
    threads: 16
    shell:
        # note: diamond messages go to stderr
        # in the output fmt, cols 6 and 7 need to be bitscore and taxid
        # for the scripts subset_blast.py and tally_abundant_hosts.py
        """
        diamond blastx \
            -d {params.diamond_nr_db} \
            -q {input} \
            -o {output} \
            -p {threads} \
            --evalue {params.diamond_nr_evalue} \
            -f 6 \
                qseqid \
                sseqid \
				pident \
                length \
                evalue \
                bitscore \
                staxids \
                stitle \
            2> {log}
        """

rule subset_diamond:
    message:
        """
        ** annotation **
        Retieving the 'best' hits for each {wildcards.sample} contig
        using the maximum bitscore
        """
    input:
        config["sub_dirs"]["annotation_dir"] + "/diamond/{sample}.diamond_blastx"
    output:
        config["sub_dirs"]["annotation_dir"] + "/diamond/{sample}_diamond_blastx.best_hits"
    shell:
        """
        {config[program_dir]}/scripts/subset_blast.py \
            -b {input} \
            -o {output}
        """

# if the host_depletion module was run, want to add those reads
# to the overall abundance tables here
# however, I don't want the host abundance file as a required input
# because then all the host_depletion rules will run even If not required
def get_host_reads(wildcards):
    if config["host_depletion"]:
        return(config["sub_dirs"]["depletion_dir"] + "/host/{}_host.blastn.abundance".format(wildcards.sample))
    else:
        return(False)


rule tally_diamond_organisms:
    message:
        """
        ** annotation **
        Calculating the most abundant species in the diamond results for each sample
        """
    input:
        blast = config["sub_dirs"]["annotation_dir"] + "/diamond/{sample}_diamond_blastx.best_hits",
        stats = config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts_subset.sorted.idxstats",
        depth = config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts_subset.sorted.depth",
        mapping = "logs/rRNA_mapping_summary.tsv",
    output:
        # producing both wide and long format tables here
        # the wide will be used for the report, and the long for plotting in ggplot
        config["sub_dirs"]["annotation_dir"] + "/diamond/{sample}_diamond_blastx.abundance",
    params:
        # adding host-specific stats here to append onto the tax results if available
        host = get_host_reads
    shell:
        """
        {config[program_dir]}/scripts/tally_organism_abundance.py \
            -b {input.blast} \
            -i {input.stats} \
            -d {input.depth} \
            --host {params.host} \
            -m {input.mapping} \
            -o {output} \
        """

# making this rule a 'checkpoint'
# because we don't know if all supertaxa will be detected in every sample
# this ensures that the DAG is re-evaluated depending on the files produced here
checkpoint sort_combine_abundances:
    message:
        """
        ** annotation **
        Sorting abundances into supertaxa and combining samples
        """
    input:
        expand(config["sub_dirs"]["annotation_dir"] + "/diamond/{sample}_diamond_blastx.abundance", sample=config["samples"])
    output:
        combined = config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance.all",
    params:
        stem = config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance_taxa",
    shell:
        # outputs a file for each kingdom present
        # e.g. if no viruses detected, the .vir file will not be created
        # also outputs a combined file with everything
        """
        {config[program_dir]}/scripts/sort_combine_abundances.py \
            -a {input} \
            -s {params.stem} \
            -o {output.combined}
        """

# this rule is like the 'intermediate' rule in the checkpoint examples
# it uses the kingdom wildcards to produce only plots for kingdoms that were detected

rule plot_abundances:
    input:
        config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance_taxa.{kingdom}",
    output:
        tsv = config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance_top10.{kingdom}.tsv",
        pdf = config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance_top10.{kingdom}.pdf",
        png = config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance_top10.{kingdom}.png",
    shell:
        """
        Rscript {config[program_dir]}/scripts/plot_tax_abundances.R \
            {input} {output.tsv} {output.pdf} {output.png}
        """

# can now use the get method to grab just the files produced

def aggregate_input(wildcards):
    # I think this just ensures that this function depends on the output of sort_combine_abundances
    # I don't actually use 'checkpoint_output' in this case
    checkpoint_output = checkpoints.sort_combine_abundances.get(**wildcards).output[0]
    # the wildcard_glob gets wildcards depending on the file names
    # thus, the 'kingdom' wildcard will only contain kingdoms that were output during
    # the sort_combine_abundances checkpoint
    return expand(config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance_top10.{kingdom}.png",
                  kingdom=glob_wildcards(config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance_taxa.{kingdom}").kingdom)


# this rule is required to make the whole thing work
# it takes the aggregate_input function, which expands the files required based
# on the kingkom wildcard (determined by looked at what files are present)
# because this rule requires these *top10.{kingdom}.png files, snakemake
# knows to run plot_abundances on whatever kingdoms are in the wildcard

rule aggregated:
    input:
        aggregate_input
    output:
        # this is really a dummy file to make the rules run
        # although, maybe I'll check through this file to determine
        # what png files to include in the report
        config["sub_dirs"]["annotation_dir"] + "/diamond/png_file_names.txt",
    shell:
        """
        echo {input} > {output}
        """

rule plot_overall_results:
    input:
        trim = "logs/trimmomatic_PE/trim_logs.summary",
        rRNA = "logs/rRNA_mapping_summary.tsv",
        combined = config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_abundance.all",
    output:
        pdf = "logs/overall_results.pdf",
        png = "logs/overall_results.png",
    shell:
        """
        Rscript {config[program_dir]}/scripts/plot_overall_results.R \
            {input.trim} \
            {input.rRNA} \
            {input.combined} \
            {output.pdf} {output.png}
        """

# want to create a ReST table for each individual sample
# that contains links to contigs per taxid

rule sort_contigs_taxid:
    input:
        contigs = config["sub_dirs"]["assembly_dir"] + "/spades/{sample}_assembly/transcripts.fasta",
        best_hits = config["sub_dirs"]["annotation_dir"] + "/diamond/{sample}_diamond_blastx.best_hits"
    output:
        directory(config["sub_dirs"]["annotation_dir"] + "/diamond/{sample}_contigs_taxid/")
    shell:
        """
        {config[program_dir]}/scripts/sort_contigs_taxid.py \
            -c {input.contigs} \
            -b {input.best_hits} \
            -o {output}
        """

rule create_abund_ReST_table:
    input:
        abund = config["sub_dirs"]["annotation_dir"] + "/diamond/{sample}_diamond_blastx.abundance",
        contig_dir = config["sub_dirs"]["annotation_dir"] + "/diamond/{sample}_contigs_taxid/"
    output:
        config["sub_dirs"]["annotation_dir"] + "/diamond/{sample}_diamond_blastx.abundance.ReST",
    shell:
        """
        {config[program_dir]}/scripts/abundance_ReST.py \
            -a {input.abund} \
            -d {input.contig_dir} \
            -o {output}
        """










