#!/usr/bin/env python

"""
Summarise log files from trimmomatic for plotting
This will summarise all log files into one
"""

import sys, os
import argparse

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("summarise trimmomatic logs")

parser.add_argument('-i', '--trimmomatic_input', type = str,
                    nargs = "*", help = "trimmomatic log files")
parser.add_argument('-o', '--summary_output', type = str,
                    help = "combined summary of trim logs")

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

output = open(args.summary_output, "w")
output.write("sample\tinput_pairs\tboth_surviving\tforward_only\treverse_only\tdropped\n")

for log_fl in args.trimmomatic_input:
    with open(log_fl) as fl:
        for line in fl:
            line = line.strip()
            if line.startswith("Input"):
                name = os.path.basename(log_fl).replace(".log", "")
                cols = line.split(":")
                input_pairs = cols[1].strip().split(" ")[0]
                both_surviving = cols[2].strip().split(" ")[0]
                F_surviving = cols[3].strip().split(" ")[0]
                R_surviving = cols[4].strip().split(" ")[0]
                dropped = cols[5].strip().split(" ")[0]
                output.write("\t".join([name, input_pairs, both_surviving, F_surviving, R_surviving, dropped]) + "\n")


