#!/usr/bin/env python

"""
Take *abundance files from tally_organism_abundance.py
and split into supertaxa file (Eukaryotes, Bacteria and Viruses).
Only output a file if the supertaxa is detected.
Also conbine the abundances from different samples into a
single file for plotting in ggplot
"""

import sys, os
import argparse

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("sort and combine abundance files")

parser.add_argument('-a', '--abundance', type = str,
                    nargs = "*", help = "abundance files from tally_organism_abundance.py")
parser.add_argument('-o', '--output', type = str,
                    help = "name for the combined output")
parser.add_argument('-s', '--stem', type = str,
                    help = "stem file name for the split outputs. This will be"
                           "appended with euk, bac or vir, depending"
                           " on the organisms detected")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.abundance is None or \
        args.stem is None:
    print("\n** a required input is missing\n"
          "** at least 1 abundance file and output and stem name is required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# create a list of results to write for each kingdom

euk_list = []
bac_list = []
vir_list = []

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
                euk_list.append(sample + "\t" + line + "\n")
            if kingdom == "Bacteria":
                bac_list.append(sample + "\t" + line + "\n")
            if kingdom == "Viruses":
                vir_list.append(sample + "\t" + line + "\n")

# create a combined output as well
output = open(args.output, "w")
output.write("\t".join(["Sample", "Taxid", "Kingdom", "Family", "Species", "Reads_Mapped", \
                            "Reads_Mapped_percent"]) + "\n")

# zip goes through two lists in parallel
for results, name in zip([euk_list, bac_list, vir_list], [".euk", ".bac", ".vir"]):
    # only write if something in the list
    if results:
        fl = open(args.stem + name, "w")
        fl.write("\t".join(["Sample", "Taxid", "Kingdom", "Family", "Species", "Reads_Mapped", \
                            "Reads_Mapped_percent"]) + "\n")
        for result in results:
            fl.write(result)
            output.write(result)






