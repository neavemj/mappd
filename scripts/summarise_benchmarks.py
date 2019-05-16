#!/usr/bin/env python

"""
Summarise all benchmark files in the snakemake benchmark directory
"""

import sys, os
import argparse

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("summarise benchmark files")

parser.add_argument('-b', '--benchmark_directory', type = str,
                    help = "benchmark directory containing snakemake benchmark files")
parser.add_argument('-o', '--output', type = str,
                    help = "combined summary of benchmark files")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.benchmark_directory is None or \
        args.output is None:
    print("\n** a required input is missing\n"
          "** the benchmark directory and output file name is required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# function to read a snakemake benchmark file
# skips the first header line, then returns a list of values from the second line
# fairly sure the benchmark files always only contain 2 lines

def read_benchmark(fl):
    with open(fl) as f:
        next(f)
        for line in f:
            line = line.strip()
            cols = line.split("\t")
            return(cols)


# walk through each directory in the benchmarks folder
# these should look like "01_trimmomatic", "02_spades", etc.
# the output gives the current directory, then each subdirectory, like this:
# ('benchmarks/', ['01_trimmomatic', '02_spades'], [])
# ('benchmarks/01_trimmomatic', [], ['prawn1.txt', 'prawn3.txt', 'prawn2.txt'])
# ('benchmarks/02_spades', [], ['prawn3.txt', 'prawn1.txt', 'prawn2.txt'])

output = open(args.output, "w")
output.write("\t".join(["module", "process", "sample", "s", "h:m:s", "max_rss", "max_vms", "max_uss", "max_pss", "io_in",
                        "io_out", "mean_load"]) + "\n")

for walk in os.walk(args.benchmark_directory):
    # will ignore the current directory (with no files)
    if walk[0] == args.benchmark_directory:
        continue

    # ok, now process each benchmark file
    for bench_fl in walk[2]:
        sample = bench_fl.rstrip(".txt")
        module = os.path.dirname(walk[0]).lstrip("./")
        process = os.path.basename(walk[0])
        file_path = os.path.join(walk[0], bench_fl)
        value_list = read_benchmark(file_path)
        output.write(module + "\t" + process + "\t" + sample + "\t" + "\t".join(value_list) + "\n")

