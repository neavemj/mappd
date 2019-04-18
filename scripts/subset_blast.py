#!/usr/bin/env python

"""
Take blast result and subset to 'best' match per contig
-outfmt '6 qseqid sseqid pident length evalue bitscore staxid salltitles'
"""

import sys
import argparse

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("subset blast results")

parser.add_argument('-b', '--blast', type = str,
                    help = "blast outfmt 6 with particular items (see above)")
parser.add_argument('-o', '--output', type = str,
                    help = "best blast hit per contig / query")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.blast is None or \
        args.output is None:
    print("\n** a required input is missing\n"
          "** a blast file and output name are required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# read through blast results and gather required info

blast_dict = {}

with open(args.blast) as fl:
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        contig = cols[0]
        bitscore = float(cols[5])
        if contig in blast_dict:
            if blast_dict[contig]["bitscore"] < bitscore:
                blast_dict[contig] = {"bitscore": bitscore, "line": line}
        else:
            blast_dict[contig] = {"bitscore": bitscore, "line": line}

# now write out these 'best' hits

output = open(args.output, "w")

for contig in blast_dict:
    output.write(blast_dict[contig]["line"] + "\n")



