#!/usr/bin/env python

"""
Take idxstats and depth files from bam and blast / diamond results and species abundance
Only works with single samples - won't combine sample files
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
import ete3_functions # used to get taxonomy info for taxids

# use argparse to grab command line arguments

parser = argparse.ArgumentParser("tally abundant hosts")

parser.add_argument('-b', '--blast', type = str,
                    help = "best blast hists with outfmt 6 with particular items (see above)")
parser.add_argument('-i', '--idxstats', type = str,
                    help = "samtools idxstats file from reads mapped back to assembly")
parser.add_argument('-d', '--depth', type = str,
                    help = "samtools depth file from reads mapped back to assembly")
parser.add_argument('-o', '--output', type = str,
                    help = "table with taxid, spp, and abundance for each sample for report")

# if no args given, print help and exit

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# check that the required arguments are provided

if args.blast is None or \
        args.idxstats is None or \
        args.depth is None or \
        args.output is None:
    print("\n** a required input is missing\n"
          "** a blast file, idxstats file, depth file and output name are required\n")
    parser.print_help(sys.stderr)
    sys.exit(1)

# create a dictionary of reads mapped from the idxstats file

idxstats_dict = {}
total_mapped_reads = 0

with open(args.idxstats) as fl:
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        contig = cols[0]
        length = cols[1]
        mapped = int(cols[2])
        total_mapped_reads += mapped
        idxstats_dict[contig] = mapped

# create a dictionary of base depth from the depth file

depth_dict = {}
total_bases = 0

with open(args.depth) as fl:
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        contig = cols[0]
        depth = int(cols[2])
        total_bases += depth
        if contig in depth_dict:
            depth_dict[contig] += depth
        else:
            depth_dict[contig] = depth

# create a dictionary for the blast results
# also keep a dictionary to seperately keep track of abundance
# will use this to sort taxids later

blast_dict = {}
abund_dict = {}

with open(args.blast) as fl:
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        contig = cols[0]
        # diamond sometimes returns more than 1 taxid
        # split by a semi-colon
        taxids = cols[6].split(";")
        for taxid in taxids:
            # for some reason taxid is occassionally blank
            # and it still gets added as a dict key!
            if taxid == "": continue
            if taxid in blast_dict:
                blast_dict[taxid]["mapped"] += idxstats_dict[contig]
                blast_dict[taxid]["bases"] += depth_dict[contig]
                abund_dict[taxid] += idxstats_dict[contig]
            else:
                blast_dict[taxid] = {"mapped": idxstats_dict[contig],
                                     "bases": depth_dict[contig]}
                abund_dict[taxid] = idxstats_dict[contig]

# will now use the ete3 functions to get taxonomy info for the taxids

tax_dict = {}

for taxid in blast_dict:
    # very occassionally a taxid is not found by ete3
    # this try-except will make it fail nicely
    try:
        superkingdom = ete3_functions.get_desired_rank(taxid, "superkingdom")
        family = ete3_functions.get_desired_rank(taxid, "family")
        species = ete3_functions.get_desired_rank(taxid, "species")
    except:
        print("taxid {} not found".format(taxid))
        continue
    tax_dict[taxid] = {"superkingdom": superkingdom, "family": family, "species": species}

# sort taxids by their abundance (mapped reads) for writing out
sorted_taxids = sorted([(value, key) for (key,value) in abund_dict.items()], reverse=True)

# now write results

output = open(args.output, "w")
output.write("\t".join(["Taxid", "Kingdom", "Family", "Species", "Reads_Mapped", "Reads_Mapped_percent", "Bases_Covered", "Bases_Covered_percent"]) + "\n")

for taxid in sorted_taxids:
    taxid = taxid[1]
    mapped_perc = (blast_dict[taxid]["mapped"] / total_mapped_reads) * 100
    bases_perc = (blast_dict[taxid]["bases"] / total_bases) * 100
    output.write("\t".join([taxid, tax_dict[taxid]["superkingdom"], tax_dict[taxid]["family"], \
    tax_dict[taxid]["species"], str(blast_dict[taxid]["mapped"]), str(mapped_perc), str(blast_dict[taxid]["bases"]), \
    str(bases_perc)]) + "\n")


