#!/usr/bin/env python

"""
Take blastx file from diamond results and add species, family and kingdom
The blast / diamond format should look like this:
            -outfmt '6 \
                qseqid \
                sseqid \
                pident \
                length \
                evalue \
                bitscore \
                staxid \
                stitle'
"""

import sys, os
import argparse
import pandas as pd
import numpy as np
import ete3_functions # used to get taxonomy info for taxids
import csv

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("add kingdoms to blastx")

parser.add_argument('-b', '--blast', type = str,
                    help = "best blast hists with outfmt 6 with particular items (see above)")
parser.add_argument('-o', '--output', type = str,
                    help = "blast output name with annotation columns added")

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


# go through blast results and add annotation info
output_list = []

with open(args.blast) as fl:
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        contig = cols[0]
        # diamond sometimes returns more than 1 taxid
        # split by a semi-colon
        # initially I ran the below loop for each taxid
        # however, this causes a doubling up of reads
        # can only run each taxid once
        taxid = cols[6].split(";")[0]
        # sometimes taxid is empty for some reason
        if taxid == "":
            superkingdom = "<taxid not found>"
            family = "<taxid not found>"
            species = "<taxid not found>"
            continue
        # now try to get classification using ete3
        # fail nicely if not detected
        try:
            superkingdom = ete3_functions.get_desired_rank(taxid, "superkingdom")
            family = ete3_functions.get_desired_rank(taxid, "family").split(",")[0]
            species = ete3_functions.get_desired_rank(taxid, "species")
        except:
            print("taxid {} not found".format(taxid))
            superkingdom = "<taxid not found>"
            family = "<taxid not found>"
            species = "<taxid not found>"
            
        # add these results to output
        output_list.append(cols +  [superkingdom, family, species])
        
# now lists of lists to file
# [['3804e5b2-0ca1-4adb-aa07-6d56e085d70f', 'KXJ07667.1', '96.6', '89', '3.45e-50', '184', '2652724', 'KXJ07667.1 hypothetical protein AC249_AIPGENE26950 [Exaiptasia diaphana]', '43.4', 'Eukaryota', 'Aiptasiidae', 'Exaiptasia diaphana'], ['011ea375-8b6a-484e-b82f-e5f4315641e9', 'CAE7942151.1', '41.3', '143', '2.65e-10', '70.9', '1628268', 'CAE7942151.1 unnamed protein product [Symbiodinium necroappetens]', '13.7', 'Eukaryota', 'Symbiodiniaceae', 'Symbiodinium necroappetens']]
with open(args.output, "w") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerows(output_list)

 
            


