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
# removed depth at this stage - could be used to calculate 'bases covered' etc
#parser.add_argument('-d', '--depth', type = str,
#                    help = "samtools depth file from reads mapped back to assembly")
parser.add_argument('-m', '--mapping', type = str,
                    help = "mapping_summary.tsv file required to get overall read numbers")
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
        args.output is None:
    print("\n** a required input is missing\n"
          "** a host table, idxstats file and output name are required\n")
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
#
#total_bases = 0
#
#with open(args.depth) as fl:
#    for line in fl:
#        line = line.strip()
#        cols = line.split("\t")
#        contig = cols[0]
#        depth = int(cols[2])
#        total_bases += depth

# get overall mRNA read numbers from mapping_summary.tsv to calculate percentages
# tricky thing is that this script runs once for each sample,
# however, the mapping file contains results for all samples

with open(args.mapping) as fl:
    sample = os.path.basename(args.output).split("_")[0]
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        name = cols[0]
        type = cols[1]
        paired_reads = cols[2]
        if name == sample and type == "mRNA_pairs":
            overall_reads = int(paired_reads) * 2

# now write results with the host taxid info

output = open(args.output, "w")
output.write("\t".join(["Taxid", "Kingdom", "Family", "Species", "Reads_Mapped", "Reads_Mapped_percent"]) + "\n")

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

        # get percentage mapped compared to all high-quality reads in dataset
        # then round to 4 decimal places
        mapped_perc = round((total_mapped_reads / overall_reads) * 100, 4)

        output.write("\t".join([taxid, superkingdom, family, species, str(total_mapped_reads),
                                str(mapped_perc)]) + "\n")

        # just get the most abundand taxid as per other modules
        # could change if more than 1 host is downloaded in the
        # host_depletion module
        break


