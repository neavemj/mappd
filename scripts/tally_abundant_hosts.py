#!/usr/bin/env python

"""
Take most abundant blast hits from multiple samples and get most abundant species
-outfmt '6 qseqid sseqid pident length evalue bitscore staxid salltitles'
"""

import sys, os
import argparse
import pandas as pd
import ete3_functions # used to get taxonomy info for taxids

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("tally abundant hosts")

parser.add_argument('-b', '--blast', type = str,
                    nargs = "*", help = "best blast hists with outfmt 6 with particular items (see above)")
parser.add_argument('-t', '--table', type = str,
                    help = "table with taxid, spp, and hits for each sample for report")
parser.add_argument('-l', '--long', type = str,
                    help = "table with taxid, spp, and hits for each sample but in long format for plotting")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.blast is None or \
        args.table is None or \
        args.long is None:
    print("\n** a required input is missing\n"
          "** a blast file and output name are required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# create a dictionary for the blast results
# and also keep a record of the most abundant taxids overall (regardless of sample)

blast_dict = {}
sample_list = []
taxid_abund = {}

for blast_fls in args.blast:
    with open(blast_fls) as fl:
        sample = os.path.basename(blast_fls).split("_")[0]
        sample_list.append(sample)
        blast_dict[sample] = {}
        for line in fl:
            line = line.strip()
            cols = line.split("\t")
            taxid = cols[6]
            # most abundant taxids by sample
            if taxid in blast_dict[sample]:
                blast_dict[sample][taxid] += 1
            else:
                blast_dict[sample][taxid] = 1
            # record most abundant taxids overall
            if taxid in taxid_abund:
                taxid_abund[taxid] += 1
            else:
                taxid_abund[taxid] = 1

# will now use the ete3 functions to get taxonomy info for the taxids

tax_dict = {}

for sample in blast_dict:
    for taxid in blast_dict[sample]:
        # no need to recalculate if another sample already has this taxid
        if taxid in tax_dict: continue
        # very occassionally a taxid is not found by ete3
        # this try-except will make it fail nicely
        try:
            phylum = ete3_functions.get_desired_rank(taxid, "phylum")
            family = ete3_functions.get_desired_rank(taxid, "family")
            species = ete3_functions.get_desired_rank(taxid, "species")
        except:
            print("taxid {} not found".format(taxid))
            continue
        tax_dict[taxid] = {"phylum": phylum, "family": family, "species": species}

# now write results in wide format for the report

output_table = open(args.table, "w")

output_table.write("\t".join(["Taxid", "Phylum", "Family", "Species"]) + "\t" + "\t".join(sample_list) + "\n")

# sort taxids by their abundance for writing out
sorted_taxids = sorted([(value, key) for (key,value) in taxid_abund.items()], reverse=True)

# the sorting produces a list of tuples, second item is the taxid
for taxid in sorted_taxids:
    taxid = taxid[1]
    # the taxid might not be in there if it failed the ete3 function bit
    if taxid in tax_dict:
        output_table.write("\t".join([taxid, tax_dict[taxid]["phylum"], tax_dict[taxid]["family"], tax_dict[taxid]["species"]]))
        for sample in sample_list:
            output_table.write("\t" + str(blast_dict[sample][taxid]))
        output_table.write("\n")

# also write out in long format for plotting in ggplot

output_long = open(args.long, "w")

output_long.write("\t".join(["Taxid", "Phylum", "Family", "Species", "Sample", "Count"]) + "\n")

for taxid in sorted_taxids:
    taxid = taxid[1]
    if taxid in tax_dict:
        for sample in sample_list:
            output_long.write("\t".join([taxid, tax_dict[taxid]["phylum"], tax_dict[taxid]["family"], \
            tax_dict[taxid]["species"], sample, str(blast_dict[sample][taxid])]) + "\n")



