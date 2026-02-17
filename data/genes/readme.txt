Folder containing information related to genes and positions.

allGenesAndIdsHg19.txt: from UCSC table browser. See selected tables in provided file.
breastCancerCausalGenes.txt: list of known breast cancer genes, constructed by hand (from COSMIC CGC website)
CCGC.tsv: all COSMIC cancer genes, downloaded from the COSMIC Cancer Gene Census, downloaded 10-04-2018.
ensemblGenesHg19: from UCSC table browser. See selected tables in provided file.
hg19_proteinCodingGenes.bed: in-house file. Source unknown...


>>TWY
Created CCGC_hg38.tsv which uses v103 of the COSMIC Cancer Gene Census (already in hg38 format) with a few column
names changed to match what is expected by inputParser.py
Created ensembl_protein_coding_genes.bed from Ensembl list of protein-coding genes on chromosomes 1-22, X, and Y
