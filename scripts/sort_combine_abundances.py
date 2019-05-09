#!/usr/bin/env python

"""
Take *abundance files from tally_organism_abundance.py
and split into supertaxa file (Eukaryotes, Bacteria and Viruses).
Also conbine the abundances from different samples into a
single file for plotting in ggplot
"""

import sys, os
import argparse

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("sort and combine abundance files")

parser.add_argument('-a', '--abundance', type = str,
                    nargs = "*", help = "abundance files from tally_organism_abundance.py")
parser.add_argument('-e', '--euk_output', type = str,
                    help = "table with Eukaryotes")
parser.add_argument('-b', '--bac_output', type = str,
                    help = "table with Bacteria")
parser.add_argument('-v', '--vir_output', type = str,
                    help = "table with Viruses")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.abundance is None or \
        args.euk_output is None or \
        args.bac_output is None or \
        args.vir_output is None:
    print("\n** a required input is missing\n"
          "** at least 1 abundance file and output names for each kingdom are required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# open and write headers for the output files

euk = open(args.euk_output, "w")
bac = open(args.bac_output, "w")
vir = open(args.vir_output, "w")

for fl in [euk, bac, vir]:
    fl.write("\t".join(["Sample", "Taxid", "Kingdom", "Family", "Species", "Reads_Mapped", \
                        "Reads_Mapped_percent", "Bases_Covered", "Bases_Covered_percent"]) + "\n")

# now read through each abundance file
# and combine the samples but split by superkingdom

for abund_fl in args.abundance:
    sample = os.path.basename(abund_fl).split("_")[0]
    with open(abund_fl) as fl:
        header = next(fl)
        for line in fl:
            line = line.strip()
            cols = line.split("\t")
            kingdom = cols[1]
            if kingdom == "Eukaryota":
                euk.write(sample + "\t" + line + "\n")
            if kingdom == "Bacteria":
                bac.write(sample + "\t" + line + "\n")
            if kingdom == "Viruses":
                vir.write(sample + "\t" + line + "\n")








