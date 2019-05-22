#!/usr/bin/env python

"""
Take best_hit file from diamond or blast search
and associated contigs from the assembly
Create a directory with contigs sorted to taxid
"""

import sys, os
import argparse
from Bio import SeqIO

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("sort contigs into taxids")

parser.add_argument('-c', '--contigs', type = str,
                    help = "contigs that were classified in the best_hit results")
parser.add_argument('-b', '--best_hits', type = str,
                    help = "best_hit classification for the contigs")
parser.add_argument('-o', '--output', type = str,
                    help = "name for sorted contig directory which will be created")


# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.contigs is None or \
        args.best_hits is None or \
        args.output is None:
    print("\n** a required input is missing\n"
          "** a contigs file, best_hist file and output dir name is required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# could be that the contig files are large sometimes
# will use biopython to index contigs to save a bit on memory

contigs = SeqIO.index(args.contigs, "fasta")

# now separate contigs into separate taxids

records_to_write = {}

with open(args.best_hits) as fl:
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        header = cols[0]
        contig = contigs[header]
        taxids = cols[6].split(";")
        for taxid in taxids:
            if taxid in records_to_write:
                records_to_write[taxid].append(contig)
            else:
                records_to_write[taxid] = [contig]

# write contigs into separate taxid files in a new directory
# if the directory exists, write into existing dir
# think this should be ok. If taxid already exists in there,
# it will be overwritten with the new data

if not os.path.isdir(args.output):
    os.mkdir(args.output)

for taxid in records_to_write:
    taxid_fl = open(args.output + "/" + taxid + ".fasta", "w")
    SeqIO.write(records_to_write[taxid], taxid_fl, "fasta")






