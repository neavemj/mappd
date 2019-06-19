#!/usr/bin/env python

"""
Take abundance files for each sample
plus directory of taxid-sorted contigs
Produce a ReST compatable table linking taxids in the
table to their associated contigs
"""

import sys, os
import argparse
import maketable

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("create ReST table from abundance files")

parser.add_argument('-a', '--abundance', type = str,
                    help = ".abundance file from a particular sample")
parser.add_argument('-d', '--taxid_contig_dir', type = str,
                    help = "directory containing taxid-sorted contigs")
parser.add_argument('-o', '--output', type = str,
                    help = "name for ReST table")


# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.abundance is None or \
        args.taxid_contig_dir is None or \
        args.output is None:
    print("\n** a required input is missing\n"
          "** an abundance file, directory of taxid-sorted contigs and output name are required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# read abundance file and sort by Viruses first, then Bacteria then Eukaryotes
# should already be sorted by abundance overall
# set mapped to false in this original dict
# to avoid missing key error if a kingdom is absent

king_dict = {"vir": {"lst": [], "mapped": False, "contigs": True},
             "bac": {"lst": [], "mapped": False, "contigs": True},
             "euk": {"lst": [], "mapped": False, "contigs": True}}
taxid_list = []

with open(args.abundance) as fl:
    header = next(fl).strip().split("\t")
    # tidying up some of the headers
    header[0] = "Sequences"
    header[4] = "Identity (%)"
    header[5] = "Reads"
    header[6] = "Reads (%)"
    sample = os.path.basename(args.abundance).split("_")[0]
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        taxid = cols[0]
        taxid_list.append(cols[0])
        # if the pident is 'mapped' means that no pident is available
        # check that here so I can add a note to the table
        pident = cols[4]
        if pident == "mapped^":
            mapped = True
        else:
            mapped = False
        # sometimes a taxid doesn't have any assembled sequence
        # e.g. if it's the host or identified using single read blasting
        if os.path.isfile("{}/{}.fasta".format(args.taxid_contig_dir, taxid)):
            # put this ReST syntax around taxid to make it link to contigs
            taxid = "`" + sample + "_" + taxid + "`_"
            contigs = True
        else:
            taxid = sample + "_" + taxid + "*"
            contigs = False
        cols[0] = taxid
        kingdom = cols[1]
        if kingdom == "Viruses":
            king_dict["vir"]["lst"].append(cols)
            king_dict["vir"]["mapped"] = mapped
            king_dict["vir"]["contigs"] = contigs
        elif kingdom == "Bacteria":
            king_dict["bac"]["lst"].append(cols)
            king_dict["bac"]["mapped"] = mapped
            king_dict["bac"]["contigs"] = contigs
        elif kingdom == "Eukaryota":
            king_dict["euk"]["lst"].append(cols)
            king_dict["euk"]["mapped"] = mapped
            king_dict["euk"]["contigs"] = contigs


def create_rest(name, hdr, king_list, contigs, mapped_note):
    # if this kingdom is detected, only need to add the header line
    # and create the ReST table
    if king_list:

        new_list = [hdr] + king_list
        rest_table = maketable.make_table(new_list)
        header_str = """

*{}*
+++++++++++++++++++++++++++++++

""".format(name)
        complete_str = header_str + rest_table

        if not contigs:
            complete_str += """

* Assembled contigs not produced

"""

        if mapped_note:
            complete_str += """

^ Mapped reads do not have an exact percent identity, although it will be high (>99%)

"""

        complete_str += """

|

"""

    # if the kingdom is absent
    else:
        complete_str = """

*{}*
+++++++++++++++++++++++++++++++

NONE DETECTED.


""".format(name)

    return(complete_str)


euk_out = create_rest("Eukaryotes", header, king_dict["euk"]["lst"], king_dict["euk"]["contigs"], king_dict["euk"]["mapped"])
bac_out = create_rest("Bacteria", header, king_dict["bac"]["lst"], king_dict["bac"]["contigs"], king_dict["bac"]["mapped"])
vir_out = create_rest("Viruses", header, king_dict["vir"]["lst"], king_dict["vir"]["contigs"], king_dict["vir"]["mapped"])

taxa_out = euk_out + bac_out + vir_out

# now need to write out the contig links to the taxids

link_str = ""

for taxid in taxid_list:
    # only add link if the fasta is present
    if os.path.isfile("{}/{}.fasta".format(args.taxid_contig_dir, taxid)):
        link_str += """
.. _{}: {}/{}.fasta
    """.format(sample + "_" + taxid, args.taxid_contig_dir, taxid)


output = open(args.output, "w")
output.write(taxa_out + "\n" + link_str)


















