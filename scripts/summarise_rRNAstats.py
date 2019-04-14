#!/usr/bin/env python

"""
Add taxonomy string from SILVA to the idxstats files
"""

import sys
import argparse

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("summarise rRNA idxstats files for plotting")

parser.add_argument('-i', '--idxstats', type = str,
                    help = ".idxstats file")
parser.add_argument('-t', '--taxstring', type = str,
                    help = ".idxstats.taxstring file")
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

with open(args.taxstring) as fl:
    for line in fl:
        line = line.strip()
        acc = line.split(" ")[0].lstrip(">")
        taxstring = " ".join(line.split(" ")[1:])
        acc_to_taxstring[acc] = taxstring

# now read through taxids file and output taxonomy info

output = open(args.output, "w")
output.write("\t".join(["accession", "mapped_reads", "taxstring"]) + "\n")

with open(args.idxstats) as fl:
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        acc = cols[0]
        taxstring = acc_to_taxstring[acc]
        mapped_reads = int(cols[2])
        output.write("\t".join([acc, str(mapped_reads), taxstring]) + "\n")




