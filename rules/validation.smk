"""
These rules will only run during validation mode
They will simulate Illumina reads based on a table
specify which genomes to use and how many reads to simulate
from each genome.
Then output a table / figure displaying the sensitivity and
specificity of the pipeline on the simulated reads

"""

import pandas as pd

virus_config = pd.read_csv("virus_config.txt", sep="\t").set_index("species_name", drop=False)
vir_genomes = list(virus_config.index.values)
host_config = pd.read_csv("host_config.txt", sep="\t").set_index("species_name", drop=False)
host_genomes = list(host_config.index.values)


rule simulate_virus_reads:
    input:
    output:
        R1 = "00_simulated/{vir_genome}.virus.1.fq",
        R2 = "00_simulated/{vir_genome}.virus.2.fq",
    params:
        species = lambda wildcards: virus_config.loc[wildcards.vir_genome, 'species_name'],
        reads_req = lambda wildcards: virus_config.loc[wildcards.vir_genome, 'reads_required'],
        file_name = lambda wildcards: virus_config.loc[wildcards.vir_genome, 'file'],
        out_stem = "00_simulated/{vir_genome}.virus."
    run:
        # sometimes want to give a fold-coverage of the genome required
        # and sometimes want to just give the actual number of reads
        # this can be enabled by putting 'x' at the end of the number
        # e.g. (4x) which makes the variable a string
        if isinstance(params.reads_req, str):
            num_flag = "-f {}".format(params.reads_req.rstrip("x"))
        else:
            num_flag = "-c {}".format(params.reads_req)
        shell(
            """
            {config[art]}/art_illumina \
                -p \
                -i {params.file_name} \
                -l 150 \
                -ss MS \
                {num_flag} \
                -m 270 \
                -s 27 \
                -o {params.out_stem}
            """)

rule extract_host_reads:
    # host reads are already simulated and I'm just grabbing the number required
    # do it this way because hosts have lots of contigs, etc and it's difficult
    # to get the exact right amount of reads
    input:
    output:
        R1 = "00_simulated/{host_genome}.host.1.fq",
        R2 = "00_simulated/{host_genome}.host.2.fq",
    params:
        species = lambda wildcards: host_config.loc[wildcards.host_genome, 'species_name'],
        reads_req = lambda wildcards: host_config.loc[wildcards.host_genome, 'reads_required'],
        R1_file = lambda wildcards: host_config.loc[wildcards.host_genome, 'R1_file'],
        R2_file = lambda wildcards: host_config.loc[wildcards.host_genome, 'R2_file'],
    run:
        num_head = int(params.reads_req) * 4
        shell(
            """
            head \
                -n {num_head} \
                {params.R1_file} > {output.R1} &&
            head \
                -n {num_head} \
                {params.R2_file} > {output.R2}
            """)

rule cat_F_reads:
    input:
        virus = expand("00_simulated/{vir_genome}.virus.1.fq", vir_genome=vir_genomes),
        host = expand("00_simulated/{host_genome}.host.1.fq", host_genome=host_genomes)
    output:
        "00_simulated/combined_sim.R1.fq"
    shell:
        """
        cat {input} > {output}
        """

rule cat_R_reads:
    input:
        virus = expand("00_simulated/{vir_genome}.virus.2.fq", vir_genome=vir_genomes),
        host = expand("00_simulated/{host_genome}.host.2.fq", host_genome=host_genomes)
    output:
        "00_simulated/combined_sim.R2.fq"
    shell:
        """
        cat {input} > {output}
        """

rule shuffle_reads:
    input:
        R1 = "00_simulated/combined_sim.R1.fq",
        R2 = "00_simulated/combined_sim.R2.fq"
    output:
        R1 = "00_simulated/combined_sim_shuf.R1.fq",
        R2 = "00_simulated/combined_sim_shuf.R2.fq"
    shell:
        """
        shuffle.sh \
            in={input.R1} \
            in2={input.R2} \
            out={output.R1} \
            out2={output.R2}
        """



# then put these cat'd reads as the sample input in my config.yaml
# this makes snakemake run these steps! cool.
# now add a rule to generate the output I want
# probably take the abundance files and generate sensitivity / specificity plots




