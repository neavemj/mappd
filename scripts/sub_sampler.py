#!/usr/bin/env python

"""
subsample a fastq file to reduce assembly time
"""

import sys
import argparse
import random
from Bio.SeqIO.QualityIO import FastqGeneralIterator    # requires Biopython

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("randomly subsample fastq files")

parser.add_argument('-1', '--forward_reads', type = str,
        nargs = "?", help = "fastq file containing forward R1 reads")
parser.add_argument('-2', '--reverse_reads', type = str,
        nargs = "?", help = "fastq file containing reverse R2 reads (leave blank if only single-end reads)")
parser.add_argument('-n', '--nReads', type = int,
        nargs = "?", default=100000, help = "number of reads required")
parser.add_argument('-o', '--output_prefix', type = str,
                    nargs = "?", help = "output prefix for the subsetted file")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.forward_reads is None or \
    args.nReads is None or \
    args.output_prefix is None:
        print("\n** a required input is missing\n"
              "** a reads file, number of reads to subset, and output prefix is required\n")
        parser.print_help(sys.stderr)
        sys.exit(1)

# figure out how many total reads are in the file

print("Scanning fastq file and calculating number of reads")

total_reads = 0

with open(args.forward_reads) as f:
    for title, seq, qual in FastqGeneralIterator(f):
        total_reads += 1

print("detected {} reads in the forward file".format(total_reads))

# create random numbers of reads to sample

reads_to_sample = set(random.sample(range(total_reads + 1), args.nReads))

print(len(reads_to_sample))

# go through forward reads again and write sub-sampled reads to file

record_number = 0
forward_ids = set()
output_forward = open(args.output_prefix + ".subset.R1.fastq", "w")

with open(args.forward_reads) as f:
    for title, seq, qual in FastqGeneralIterator(f):
        record_number += 1
        if record_number in reads_to_sample:
            forward_ids.add(title)
            output_forward.write("@{}\n{}\n+\n{}\n".format(title, seq, qual))




#reverse_handle = open(args.output_prefix + ".R2.fq", "w")
