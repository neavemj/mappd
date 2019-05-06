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
        Using Diamond blastx to compare contigs to the nr database
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
        "benchmarks/diamond/{sample}.txt"
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
        Retieving the 'best' hits for each {wildcards.sample} contig
        using the maximum bitscore
        """
    input:
        config["sub_dirs"]["annotation_dir"] + "/diamond/{sample}.diamond_blastx"
    output:
        config["sub_dirs"]["annotation_dir"] + "/diamond/{sample}.diamond_blastx.best_hits"
    shell:
        """
        {config[program_dir]}/scripts/subset_blast.py \
            -b {input} \
            -o {output}
        """

# could use idxstats to get reads mapped to each SPAdes contig
# then could combine this with the taxonomy results
# to output a table showing abundant organisms with how many reads mapped

rule tally_diamond_organisms:
    message:
        """
        Calculating the most abundant species in the diamond results for all samples
        """
    input:
        config["sub_dirs"]["annotation_dir"] + "/diamond/{sample}.diamond_blastx.best_hits"
    output:
        # producing both wide and long format tables here
        # the wide will be used for the report, and the long for plotting in ggplot
        wide = config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_tax.wide",
        long = config["sub_dirs"]["annotation_dir"] + "/diamond/diamond_blastx_tax.long",
    shell:
        """
        {config[program_dir]}/scripts/tally_blast_organisms.py \
            -b {input} \
            -t {output.wide} \
            -l {output.long}
        """



