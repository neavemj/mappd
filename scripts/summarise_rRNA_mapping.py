#!/usr/bin/env python

"""
Summarise bowtie results for the rRNA databases
This will summarise multiple sample results into a single file
"""

import sys, os
import argparse

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("Summarise bowtie rRNA mapping results")

parser.add_argument('-l', '--lsu', type = str,
                    nargs = "*", help = "LSU log files")
parser.add_argument('-s', '--ssu', type = str,
                    nargs = "*", help = "SSU log files")
parser.add_argument('-o', '--output', type = str,
                    help = "summary of rRNA mapping")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.lsu is None or \
        args.ssu is None:
    print("\n** a required input is missing\n"
          "** both rRNA (SSU and LSU) mapping log files are required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# read through each log file and extract required info

def summarise_bowtie(log_fl):
    sample = os.path.basename(log_fl).split(".")[0]
    line_count = 0
    with open(log_fl) as fl:
        for line in fl:
            line_count += 1
            line = line.strip()
            if line_count == 1:
                total_reads = line.split()[0]
            if "pairs aligned 0 times concordantly or discordantly" in line:
                unaligned = line.split()[0]
    return(sample, total_reads, unaligned)


sample_to_mapping = {}

for log_fl in args.lsu:
    results = summarise_bowtie(log_fl)
    sample_to_mapping[results[0]] = {"LSU": {"total": results[1], "unaligned": results[2]}}

for log_fl in args.ssu:
    results = summarise_bowtie(log_fl)
    sample_to_mapping[results[0]]["SSU"] = {"total": results[1], "unaligned": results[2]}

# looks good now write out file for plotting
# have to calculate the actual aligned reads weirdly like this due to bowtie output
# because I'm removing mapped pairs, including if only 1 pair mapped,
# it doesn't appear in the output

output = open(args.output, "w")
output.write("\t".join(["Sample", "Type", "Paired_Reads"]) + "\n")

for sample in sample_to_mapping:
    LSU_total = int(sample_to_mapping[sample]["LSU"]["total"])
    SSU_total = int(sample_to_mapping[sample]["SSU"]["total"])

    LSU_pairs = LSU_total - SSU_total

    SSU_pairs = SSU_total - int(sample_to_mapping[sample]["SSU"]["unaligned"])

    mRNA_pairs = sample_to_mapping[sample]["SSU"]["unaligned"]

    output.write("\t".join([sample, "rRNA_LSU", str(LSU_pairs)]) + "\n")
    output.write("\t".join([sample, "rRNA_SSU", str(SSU_pairs)]) + "\n")
    output.write("\t".join([sample, "mRNA_pairs", str(mRNA_pairs)]) + "\n")











