# need this to ensure the SAMPLES can be accessed anywhere

from os.path import join, abspath
from snakemake.io import glob_wildcards

def get_samples(raw_dir, illumina_machine):
    raw_dir = abspath(raw_dir)

    if illumina_machine == "nextseq":
        pattern = join(raw_dir, "{sample}_S{run}_R1_001.fastq.gz")
    elif illumina_machine == "miseq":
        pattern = join(raw_dir, "{sample}_S{run}_L001_R1_001.fastq.gz")
    else:
        raise ValueError(f"Unknown illumina_machine: {illumina_machine}")

    samples = glob_wildcards(pattern).sample
    #print(f"DEBUG: Matched samples: {samples}")
    return samples
