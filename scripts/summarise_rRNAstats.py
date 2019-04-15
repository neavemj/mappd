#!/usr/bin/env python

"""
Add taxonomy string from SILVA to the idxstats files
This will summarise sample results into a single file
"""

import sys, os
import argparse

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("summarise rRNA idxstats files for plotting")

parser.add_argument('-i', '--idxstats', type = str,
                    nargs = "*", help = ".idxstats files")
parser.add_argument('-t', '--taxstring', type = str,
                    nargs = "*", help = ".idxstats.taxstring files")
parser.add_argument('-o', '--output', type = str,
                    help = "summary of rRNA mapping")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.idxstats is None or \
        args.taxstring is None:
    print("\n** a required input is missing\n"
          "** both the idxstats and taxstring files are required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# create a lookup dictionary for the acc_to_taxstring conversion

acc_to_taxstring = {}

for taxstring_fls in args.taxstring:
    with open(taxstring_fls) as fl:
        for line in fl:
            line = line.strip()
            acc = line.split(" ")[0].lstrip(">")
            taxstring = " ".join(line.split(" ")[1:])
            acc_to_taxstring[acc] = taxstring

# now read through idxstats files and output taxonomy info

output = open(args.output, "w")
output.write("\t".join(["Sample", "Accession", "rRNA_Type", "Mapped_Reads", "Taxonomy_String", "Species"]) + "\n")

for idxstats_fls in args.idxstats:
    with open(idxstats_fls) as fl:
        for line in fl:
            line = line.strip()
            cols = line.split("\t")
            acc = cols[0]
            taxstring = acc_to_taxstring[acc]
            species = taxstring.split(";")[-1]
            mapped_reads = int(cols[2])
            sample = os.path.basename(idxstats_fls).split("_")[0]
            rRNA = os.path.basename(idxstats_fls).split("_")[1].split(".")[0]
            output.write("\t".join([sample, acc, rRNA, str(mapped_reads), taxstring, species]) + "\n")




