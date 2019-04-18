#!/usr/bin/env python

"""
Gather a certain number of contigs larger than a given threshold
"""
import sys
import argparse
from Bio import SeqIO # biopython required

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("gather contigs")

parser.add_argument('-c', '--contigs', type = str,
                    help = "fasta file containing contigs")
parser.add_argument('-s', '--size', type = str,
                    help = "minimum contig length to output")
parser.add_argument('-n', '--num', type = str,
                    help = "maximum number of contigs to output"
                           "can be set to a very large number to output all")
parser.add_argument('-o', '--output', type = str,
                    help = "output of subset of contigs")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.contigs is None or \
        args.size is None or \
        args.num is None or \
        args.output is None:
    print("\n** a required input is missing\n"
          "** a contig file, min size, max number and output name are required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# read through contigs and output as required

seqs_to_write = []
count = 1

with open(args.contigs) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if len(record.seq) > int(args.size) and count <= int(args.num):
            count += 1
            seqs_to_write.append(record)

SeqIO.write(seqs_to_write, open(args.output, "w"), "fasta")












