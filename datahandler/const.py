"""Human net is a probabilistic functional gene network of 18,714 validated protein-encoding 
genes of Homo sapiens (by NCBI March 2007), constructed by a modified Bayesian integration
of 21 types of 'omics' data from multiple organisms, with each data type weighted according
to how well it links genes that are known to function together in H. sapiens. Each interaction
in HumanNet has an associated log-likelihood score (LLS) that measures the probability
of an interaction representing a true functional linkage between two genes.

v.1: # genes = 16,243;	# linkages = 476,399;   coverage of 18,714 validated protein-coding loci = 86.8%

For more check:
Prioritizing candidate disease genes by network-based boosting of genome-wide association data
Insuk Lee, U. Martin Blom, Peggy I. Wang, Jung Eun Shin, and Edward M. Marcotte
Genome Research 21(7):1109-21 (2011)
"""
GENE_GENE_SIMILARITY_URL: str = (
    "http://www.functionalnet.org/humannet/HumanNet.v1.join.txt"
)
MESH_URL: str = "https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/"
OMIM_URL: str = "https://data.omim.org/downloads/"
API_KEY: str = "OMIM_API_KEY"
