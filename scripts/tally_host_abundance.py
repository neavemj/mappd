#!/usr/bin/env python

"""
Take idxstats and depth files from host mapping and summarise abundance
Only works with single samples - won't combine sample files
The host table should be in wide format
"""

import sys, os
import argparse
import pandas as pd
import ete3_functions # used to get taxonomy info for taxids

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("summarise host abundance")

parser.add_argument('-w', '--wide', type = str,
                    help = "most abundant host table in wide format")
parser.add_argument('-i', '--idxstats', type = str,
                    help = "samtools idxstats file from reads mapped back to assembly")
parser.add_argument('-d', '--depth', type = str,
                    help = "samtools depth file from reads mapped back to assembly")
parser.add_argument('-o', '--output', type = str,
                    help = "table with taxid, spp, and abundance for each sample for report")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.wide is None or \
        args.idxstats is None or \
        args.depth is None or \
        args.output is None:
    print("\n** a required input is missing\n"
          "** a host table, idxstats file, depth file and output name are required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# get reads mapped from the idxstats file

total_mapped_reads = 0

with open(args.idxstats) as fl:
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        contig = cols[0]
        length = cols[1]
        mapped = int(cols[2])
        total_mapped_reads += mapped

# get base depth from the depth file

total_bases = 0

with open(args.depth) as fl:
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        contig = cols[0]
        depth = int(cols[2])
        total_bases += depth

# now write results with the host taxid info

output = open(args.output, "w")
output.write("\t".join(["Taxid", "Kingdom", "Family", "Species", "Reads_Mapped", "Reads_Mapped_percent", "Bases_Covered", "Bases_Covered_percent"]) + "\n")

with open(args.wide) as fl:
    next(fl)
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        taxid = cols[0]
        # will now use the ete3 functions to get taxonomy info for the taxids
        try:
            superkingdom = ete3_functions.get_desired_rank(taxid, "superkingdom")
            family = ete3_functions.get_desired_rank(taxid, "family").split(",")[0]
            species = ete3_functions.get_desired_rank(taxid, "species")
        except:
            print("taxid {} not found".format(taxid))
            continue

        # need to calculate this better
        # currently (in tally_organism_abundance.py) it's calculated based on
        # total reads mapped to the assembly
        # maybe should be bases on total reads in the entire dataset?
        # could then incorporate that figure here
        mapped_perc = "100"
        bases_perc = "100"

        output.write("\t".join([taxid, superkingdom, family, species, str(total_mapped_reads),
                                mapped_perc, str(total_bases), bases_perc]) + "\n")

        # just get the most abundand taxid as per other modules
        # could change if more than 1 host is downloaded in the
        # host_depletion module
        break


