#!/usr/bin/env python

"""
Take idxstats and depth files from bam and blast / diamond results and species abundance
Only works with single samples - won't combine sample files
Can take a pre-calculated host abundance file to add to the results calculated here
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

parser = argparse.ArgumentParser("tally abundant organisms")

parser.add_argument('-b', '--blast', type = str,
                    help = "best blast hists with outfmt 6 with particular items (see above)")
parser.add_argument('-i', '--idxstats', type = str,
                    help = "samtools idxstats file from reads mapped back to assembly")
parser.add_argument('-d', '--depth', type = str,
                    help = "samtools depth file from reads mapped back to assembly")
parser.add_argument('--host', type = str,
                    help = "a pre-calculated host abundance file to add to these stats. "
                            "Can leave blank if not required.")
parser.add_argument('-m', '--mapping', type = str,
                    help = "mapping_summary.tsv file required to get overall read numbers")
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
        args.mapping is None or \
        args.output is None:
    print("\n** a required input is missing\n"
          "** a blast file, idxstats file, depth file, mapping file and output name are required\n")
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
        # initially I ran the below loop for each taxid
        # however, this causes a doubling up of reads
        # can only run each taxid once
        taxid = cols[6].split(";")[0]

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
        family = ete3_functions.get_desired_rank(taxid, "family").split(",")[0]
        species = ete3_functions.get_desired_rank(taxid, "species")
    except:
        print("taxid {} not found".format(taxid))
        continue
    tax_dict[taxid] = {"superkingdom": superkingdom, "family": family, "species": species}

# sort taxids by their abundance (mapped reads) for writing out
sorted_taxids = sorted([(value, key) for (key,value) in abund_dict.items()], reverse=True)

# get overall mRNA read numbers from mapping_summary.tsv to calculate percentages
# tricky thing is that this script runs once for each sample,
# however, the mapping file contains results for all samples

with open(args.mapping) as fl:
    sample = os.path.basename(args.output).split("_")[0]
    for line in fl:
        line = line.strip()
        cols = line.split("\t")
        name = cols[0]
        type = cols[1]
        paired_reads = cols[2]
        if name == sample and type == "mRNA_pairs":
            overall_reads = int(paired_reads) * 2

if not overall_reads:
    print("could not determine overall reads from mapping_summary.tsv file")
    print("outputting 0")
    overall_reads = 0

# now write results
output = open(args.output, "w")
output.write("\t".join(["Taxid", "Kingdom", "Family", "Species", "Reads_Mapped", "Reads_Mapped_percent"]) + "\n")

# if a pre-calculated host abundance file is created, will be added first
# should be the most abundant taxid
if not args.host == "False":
    with open(args.host) as fl:
        next(fl)
        for line in fl:
            output.write(line)

# now add other taxids in order of abundance

for taxid in sorted_taxids:
    taxid = taxid[1]
    # get percentage mapped compared to all high-quality reads in dataset
    # then round to 4 decimal places
    mapped_perc = round((blast_dict[taxid]["mapped"] / overall_reads) * 100, 4)
    # note: not writing out bases covered at this stage
    bases_perc = (blast_dict[taxid]["bases"] / total_bases) * 100
    output.write("\t".join([taxid, tax_dict[taxid]["superkingdom"], tax_dict[taxid]["family"], \
    tax_dict[taxid]["species"], str(blast_dict[taxid]["mapped"]), str(mapped_perc)]) + "\n")



