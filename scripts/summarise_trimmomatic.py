#!/usr/bin/env python

"""
Summarise the log file from trimmomatic for plotting
This will only process one log file at a time; should work better with snakemake
I will combine the output files at the plotting stage
"""

import sys
import argparse

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("summarise trimmomatic log")

parser.add_argument('-i', '--trimmomatic_input', type = str,
                    nargs = "?", help = "fastq file containing forward R1 reads")
parser.add_argument('-o', '--summary_output', type = str,
                    nargs = "?", help = "fastq file containing reverse R2 reads (leave blank if only single-end reads)")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.trimmomatic_input is None or \
        args.summary_output is None:
    print("\n** a required input is missing\n"
          "** a trimmomatic log and output file name is required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# read through the file and grab the important information

with open(args.trimmomatic_input) as fl:
    output = open(args.summary_output, "w")
    output.write("input_pairs\tboth_surviving\tforward_only\treverse_only\tdropped\n")
    for line in fl:
        line = line.strip()
        if line.startswith("Input"):
            cols = line.split(":")
            input_pairs = cols[1].strip().split(" ")[0]
            both_surviving = cols[2].strip().split(" ")[0]
            F_surviving = cols[3].strip().split(" ")[0]
            R_surviving = cols[4].strip().split(" ")[0]
            dropped = cols[5].strip().split(" ")[0]

            output.write("\t".join([input_pairs, both_surviving, F_surviving, R_surviving, dropped]))


