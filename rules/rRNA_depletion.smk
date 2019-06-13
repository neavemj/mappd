"""
Ribosomal RNA depletion rules

These rules will map reads to the SILVA LSU and SSU rRNA databases
and split the reads that don't match
Another advantage of bbmap is that the fastq headers remain intact
This is important for downstream tools, such as trinity

"""

rule bbmap_to_LSU:
    message:
        """
        ** rRNA_depletion **
        Mapping cleaned {wildcards.sample} reads to the SILVA LSU rRNA database
        """
    input:
        R1 = config["sub_dirs"]["trim_dir"] + "/{sample}_1P.fastq.gz",
        R2 = config["sub_dirs"]["trim_dir"] + "/{sample}_2P.fastq.gz"
    output:
        R1 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU_depleted_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU_depleted_2P.fastq"
    params:
        silva_LSU_db = config['silva_LSU_db'],
        max_memory = config['bbmap_max_memory'],
    log:
        "logs/bbmap_LSU/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["depletion_dir"] + "/bbmap_LSU/{sample}.txt"
    threads: 16
    shell:
        """
        bbmap.sh \
            in1={input.R1} \
            in2={input.R2} \
            outu1={output.R1} \
            outu2={output.R2} \
            threads={threads} \
            {params.max_memory} \
            path={params.silva_LSU_db} \
            > {log}
        """

rule bbmap_to_SSU:
    message:
        """
        ** rRNA_depletion **
        Mapping LSU-depleted {wildcards.sample} reads to the SILVA SSU rRNA database
        """
    input:
        R1 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU_depleted_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU_depleted_2P.fastq"
    output:
        R1 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_mRNA_1P.fastq",
        R2 = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_mRNA_2P.fastq"
    params:
        silva_SSU_db = config['silva_SSU_db'],
        max_memory = config['bbmap_max_memory'],
    log:
        "logs/bbmap_SSU/{sample}.log"
    benchmark:
        "benchmarks/" + config["sub_dirs"]["depletion_dir"] + "/bbmap_SSU/{sample}.txt"
    threads: 16
    shell:
        """
        bbmap.sh \
            in1={input.R1} \
            in2={input.R2} \
            outu1={output.R1} \
            outu2={output.R2} \
            threads={threads} \
            {params.max_memory} \
            path={params.silva_SSU_db} > {log}
        """

rule summarise_sample_rRNA:
    message:
        """
        ** rRNA_depletion **
        Calculating number of reads mapped to LSU and SSU databases
        """
    input:
        trimmed = config["sub_dirs"]["trim_dir"] + "/{sample}_1P.fastq.gz",
        LSU_depleted = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_LSU_depleted_1P.fastq",
        SSU_depleted = config["sub_dirs"]["depletion_dir"] + "/rRNA/{sample}_mRNA_1P.fastq",
    output:
        "logs/rRNA_summary/{sample}_rRNA_mapping_summary.tsv"
    shell:
        # first calculates how many reads in each file (divide by 2 gives all reads)
        # then uses these numbers to calculated actual LSU and SSU values
        # then write this info to file
        # have to do it ugly like this because echo inserts a space otherwise
        """
        trimmed_reads=$(expr $(zcat {input.trimmed} | wc -l) / 2)
        LSU_depleted_reads=$(expr $(cat {input.LSU_depleted} | wc -l) / 2)
        mRNA=$(expr $(cat {input.SSU_depleted} | wc -l) / 2)

        LSU_reads=$(expr $trimmed_reads - $LSU_depleted_reads)
        SSU_reads=$(expr $LSU_depleted_reads - $mRNA)

        echo -e \
            {wildcards.sample}'\trRNA_LSU\t'${{LSU_reads}}'\n'{wildcards.sample}'\trRNA_SSU\t'${{SSU_reads}}'\n'{wildcards.sample}'\tmRNA_reads\t'${{mRNA}}'\n'\
                > {output}
        """

rule summarise_rRNA_mapping:
    input:
        expand("logs/rRNA_summary/{sample}_rRNA_mapping_summary.tsv", sample=config["samples"])
    output:
        "logs/rRNA_mapping_summary.tsv"
    shell:
        """
        cat {input} > {output}
        """












