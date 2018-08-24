
configfile: "config.yaml"

rule all:
    input:
        "test.kraken.txt"

rule count_reads:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "{sample}.read_count"
    shell:
        "grep -c '@' {input} > {output}"

# this kraken run doesn't work
# think because the database was built for kraken1 (not v2)

rule run_kraken:
    input:
        "data/test.fasta"
    output:
        "test.kraken.txt"
    shell:
        "/data/nea040/software/kraken2-2.0.7-beta/kraken2 --db /datastore/nea040/kraken_db/minikraken_20171101_4GB_dustmasked/ {input} > {output}"


