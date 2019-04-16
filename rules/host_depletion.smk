"""
Host depletion rules

These rules will map reads to host genetic data
and split the reads into those that match and those that don't
These rules rely on the results from 'rRNA_depletion.smk'
in order to identy the most likely host species
Note: these only use results from the SSU database
as these are more accurate because more people
sequence the SSU (e.g. 18S)
"""


rule associate_hostTaxid_genbank:
    message:
        """
        Retreiving genbank ids for the most abundant potential host species
        Currently this is grepping all 10 most abundant
        Maybe I should do 1 at a time - might end up with lots of sequence
        """
    input:
        config["sub_dirs"]["depletion_dir"] + "/rRNA/SSU.idxstats.summary.tsv"
    output:
        config["sub_dirs"]["depletion_dir"] + "/host/host1_nucl_gb.ids"
    params:
        acc_to_taxids = config["acc_to_taxids"],
        # need to do this to ensure only the string matches
        sed_pat = r"s/\(.*\)/\t\1\t/g"
    shell:
        # cuts the second field (taxid), removes the header,
        # adds a tab before and after the taxid to ensure a
        # clean match
        """
        grep \
            -f <(cut -f 2 {input} | tail -n+2 | sed "{params.sed_pat}") \
            {params.acc_to_taxids} \
            > {output}
        """


