# these functions use the ete3 package to manipulate the NCBI taxonomy
# this is used in the host depletion step and for summarising the annotations

# the first time this is used, it will download the NCBI taxonomy database into your home directory
import os
from ete3 import NCBITaxa

# Custom location of the ete3 taxon database
# Default locaition is the ~/.etetoolkit directory
taxadb = os.getenv('ETE3_TAXONDB')
print("=== dbname ===")
print(taxadb)
ncbi = NCBITaxa(taxadb)


# helper function to return the desired rank from a taxid


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

