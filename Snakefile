
configfile: "config.yaml"

rule all:
    input:
        "novel.read_count"

rule count_reads:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "{sample}.read_count"
    shell:
        "grep -c '@' {input} > {output}"
