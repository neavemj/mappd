#!/usr/bin/env python

"""
Take *abundance files from tally_organism_abundance.py
and associated contigs from the assembly
Create a restructuredtext table with links between species
and the contigs identified to the species
Sort results by viruses, bacteria, then eukaryotes with a
restructuredtext line between them?
"""

import sys, os
import argparse

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("create ReST table from abundance file with contig links")

parser.add_argument('-a', '--abundance', type = str,
                    help = "abundance files from tally_organism_abundance.py")
parser.add_argument('-b', '--best_hits', type = str,
                    help = "best_hits file to get contig names")
parser.add_argument('-c', '--contigs', type = str,
                    help = "contigs that were classified for the abundance files")
parser.add_argument('-o', '--output', type = str,
                    help = "name for ReST table")
parser.add_argument('-t', '--output_contigs', type = str,
                    help = "name for directory to store taxid-sorted contigs")


# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.abundance is None or \
        args.contigs is None or \
        args.output is None or \
        args.output_contigs is None:
    print("\n** a required input is missing\n"
          "** an abundance file, contigs file and output names are required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)
