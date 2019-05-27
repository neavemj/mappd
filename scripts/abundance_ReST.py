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

vir_list = []
bac_list = []
euk_list = []
taxid_list = []

with open(args.abundance) as fl:
    header = next(fl).strip().split("\t")
    # tidying up some of the headers
    header[0] = "Sequences"
    header[4] = "Reads"
    header[5] = "Reads (%)"
    sample = os.path.basename(args.abundance).split("_")[0]
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        taxid_list.append(cols[0])
        # put this ReST syntax around taxid to make it link to contigs
        taxid = "`" + sample + "_" + cols[0] + "`_"
        cols[0] = taxid
        kingdom = cols[1]
        if kingdom == "Viruses":
            vir_list.append(cols)
        elif kingdom == "Bacteria":
            bac_list.append(cols)
        elif kingdom == "Eukaryota":
            euk_list.append(cols)

# sometimes the sorting gets a bit out due to the host addition
# sort the taxa lists by the 5th element (number of reads)
vir_list = sorted(vir_list, key=lambda x: int(x[4]), reverse=True)
bac_list = sorted(bac_list, key=lambda x: int(x[4]), reverse=True)
euk_list = sorted(euk_list, key=lambda x: int(x[4]), reverse=True)


def create_rest(name, hdr, king_list):
    # if this kingdom is detected, only need to add the header line
    # and create the ReST table
    if king_list:

        new_list = [hdr] + king_list
        rest_table = maketable.make_table(new_list)
        header_str = """

*{}*
+++++++++++++++++++++++++++++++

""".format(name)
        complete_str = header_str + rest_table + """

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


euk_out = create_rest("Eukaryotes", header, euk_list)
bac_out = create_rest("Bacteria", header, bac_list)
vir_out = create_rest("Viruses", header, vir_list)

taxa_out = euk_out + bac_out + vir_out

# now need to write out the contig links to the taxids

link_str = ""

for taxid in taxid_list:
    link_str += """
.. _{}: {}/{}.fasta
""".format(sample + "_" + taxid, args.taxid_contig_dir, taxid)


output = open(args.output, "w")
output.write(taxa_out + "\n" + link_str)


















