#!/bin/bash 
#SBATCH --job-name mappd 
#SBATCH --nodes 1 
#SBATCH --ntasks-per-node 1 
#SBATCH --cpus-per-task 8 
#SBATCH --mem 32gb 
#SBATCH --time 08:00:00

snakemake -s mappd.snakefile -j 8
