#!/usr/bin/env python

"""
Summarise the idxstats files from the SSU and LSU rRNA databases
This will summarise both rRNA types into one
"""

import sys
import argparse
# will use the ete3 package to manipulate the NCBI taxonomy
# the first time this is used, it will download the NCBI taxonomy database into your home directory
from ete3 import NCBITaxa

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("summarise rRNA idxstats files for plotting")

parser.add_argument('-l', '--LSU_idxstats', type = str,
                    help = "LSU.idxstats file")
parser.add_argument('-s', '--SSU_idxstats', type = str,
                    help = "SSU.idxstats file")
parser.add_argument('-t', '--LSU_taxids', type = str,
                    help = "LSU.idxstats.taxid file")
parser.add_argument('-r', '--SSU_taxids', type = str,
                    help = "SSU.idxstats.taxid file")

parser.add_argument('-o', '--output', type = str,
                    help = "combined summary of trim logs")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.LSU_idxstats is None or \
        args.SSU_idxstats is None:
    print("\n** a required input is missing\n"
          "** both the LSU and SSU idxstats files are required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# create a lookup dictionary for the acc_to_taxids conversion

acc_to_taxids = {}

with open(args.LSU_taxids) as fl:
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        acc = cols[0]
        taxid = cols[2]
        acc_to_taxids[acc] = taxid

with open(args.SSU_taxids) as fl:
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        acc = cols[0]
        taxid = cols[2]
        acc_to_taxids[acc] = taxid

# helper function to return the desired rank from a taxid

ncbi = NCBITaxa()

def get_desired_rank(taxid, desired_rank):
    lineage = ncbi.get_lineage(taxid)
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank,taxid) for (taxid, rank) in lineage2ranks.items())
    specific_taxid = ranks2lineage.get(desired_rank, '<not present>')
    if specific_taxid != '<not present>':
        return(list(ncbi.get_taxid_translator([specific_taxid]).values())[0])
    else:
        return('<not present>')

# now read through taxids file and create dict of summarised counts
# if the counts are the same taxid, I'll add them together

LSU_counts_dict = {}

with open(args.LSU_idxstats) as fl:
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        acc = cols[0].split(".")[0]
        if acc in acc_to_taxids:
            taxid = acc_to_taxids[acc]
        else:
            print("warning couldn't find taxid for: {}. Ignoring..".format(acc))
        mapped_reads = int(cols[2])
        if taxid in LSU_counts_dict:
            LSU_counts_dict[taxid] += mapped_reads
        else:
            LSU_counts_dict[taxid] = mapped_reads

SSU_counts_dict = {}

with open(args.SSU_idxstats) as fl:
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        acc = cols[0].split(".")[0]
        if acc in acc_to_taxids:
            taxid = acc_to_taxids[acc]
        else:
            print("warning couldn't find taxid for: {}. Ignoring..".format(acc))
        mapped_reads = int(cols[2])
        if taxid in SSU_counts_dict:
            SSU_counts_dict[taxid] += mapped_reads
        else:
            SSU_counts_dict[taxid] = mapped_reads

# now write this information for plotting ggplot

output = open(args.output, "w")
output.write("\t".join(["rRNA", "taxid", "mapped_reads", "phylum", "family", "species"]) + "\n")

for taxid in LSU_counts_dict:
    try:
        phylum = get_desired_rank(taxid, "phylum")
        family = get_desired_rank(taxid, "family")
        species = get_desired_rank(taxid, "species")
    except:
        print("Couldn't find taxonomy string for taxid: {}".format(taxid))
    output.write("\t".join(["LSU", taxid, str(LSU_counts_dict[taxid]), phylum, family, species]) + "\n")

for taxid in SSU_counts_dict:
    try:
        phylum = get_desired_rank(taxid, "phylum")
        family = get_desired_rank(taxid, "family")
        species = get_desired_rank(taxid, "species")
    except:
        print("Couldn't find taxonomy string for taxid: {}".format(taxid))
    output.write("\t".join(["SSU", taxid, str(SSU_counts_dict[taxid]), phylum, family, species]) + "\n")

















