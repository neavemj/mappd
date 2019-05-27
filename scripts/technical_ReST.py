#!/usr/bin/env python

"""
Take technical files, such as start, end, software versions, etc.
and produce a ReST compatable table
"""

import sys, os
import argparse
import maketable
import json

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("create ReST table from technical files")

parser.add_argument('-c', '--config', type = str,
                    help = "config file")
parser.add_argument('-s', '--start', type = str,
                    help = "start log")
parser.add_argument('-e', '--end', type = str,
                    help = "end log")
parser.add_argument('-f', '--software', type = str,
                    help = "software list")
parser.add_argument('-o', '--output', type = str,
                    help = "name for technical ReST table")


# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.start is None or \
        args.end is None or \
        args.software is None or \
        args.config is None or \
        args.output is None:
    print("\n** a required input is missing\n"
          "** a start log, end log, software list, config file and output name are required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

## Sample Info ##

# read the sample config info into a proper dictionary
# it is initally called in as a simple string (that looks like a dict)
# apparently json.loads will do the job of converting the str to dict
# but it needs double quotes rather than single quotes

sample_dict = json.loads(args.config.replace("'", '"'))
sample_list = [["Sample", "Files"]]

# I'd like to display the files in a sample in multiple rows
# otherwise a single row gets really crowded
for sample in sample_dict:
    # the star is a dummy character that I can split by later
    # the dash makes bullet points in ReST
    files = "\t".join(["*- " + file for file in sample_dict[sample]])
    sample_list.append([sample, files])

# make the string ReST table now
sample_rest = maketable.make_table(sample_list)

# now need to hack this string to put multiple files on different rows
# first figure out the length of the first column (plus 1 for space)
col1_len = len(sample_rest.split(" ")[0]) + 1
# calculate how many spaces I need to push the second file to the correct column
col1_spaces = " " * col1_len
# use this info to replace the * with a new line plus the required number of spaces
sample_rest2 = sample_rest.replace("*", "\n" + col1_spaces)

## start and end Info ##

date_list = [["Checkpoint", "Date"]]

start_list = [start.split("_") for start in open(args.start)][0]
end_list = [end.split("_") for end in open(args.end)][0]

date_list.append(start_list)
date_list.append(end_list)

date_rest = maketable.make_table(date_list)

## software versions ##

software_list = [["Software", "Version"]]

softwares = [software.strip().split("\t") for software in open(args.software)]

software_list += softwares

software_rest = maketable.make_table(software_list)

# combine technical tables into one for writing

rest_table = date_rest + "\n\n|\n\n" + sample_rest2 + "\n\n|\n\n" + software_rest

output = open(args.output, "w")
output.write(rest_table)
















