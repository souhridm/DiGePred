The features used to characterize each gene pair to train DiGePred were:
common_pathways: Jaccard similarity of pathways associated with both genes. (https://www.genome.jp/kegg/, https://reactome.org)
common_phenotypes: Jaccard similarity of phenotypes associated with both genes. (https://hpo.jax.org/app/)
Co-expression_coefficient: mutual coexpression rank for the two genes. (https://coxpresdb.jp)
PPI_network_dist: distance on the protein-protein interaction (PPI) network. (https://genome.ucsc.edu/goldenPath/help/hgGeneGraph.html)
PWY_network_dist:	distance on the biochemical pathways interaction network. (https://genome.ucsc.edu/goldenPath/help/hgGeneGraph.html)
Txt_network_dist:	distance on the literature-mined interaction network. (https://genome.ucsc.edu/goldenPath/help/hgGeneGraph.html)
LoFintolerance: Quadratic mean of genes' sensitivity to loss-of-function variants. (https://www.deciphergenomics.org)
Haploinsufficiency: Quadratic mean of genes' resistance to haploinsufficiency. (https://www.deciphergenomics.org)
Age: Quadratic mean of gene/protein ages. (https://proteinhistorian.docpollard.org)
dN/dS: Quadratic mean of selection pressure on genes. (http://www.h-invitational.jp/evola/search.html)
Essentiality: Quadratic mean of organism's tolerance to loss of genes' funtionality. (https://v3.ogee.info/#/home)
#ofpathways_combined: Quadratic mean of number of pathways associated with each gene.
#ofPhenotypeCodes_combined: Quadratic mean of number of phenotypes associated with each gene.
#ofNeighborsPPI_combined: Quadratic mean of number of neighbors in the PPI network, associated with each gene.
#ofNeighborsPWY_combined: Quadratic mean of number of neighbors in the pathways network, associated with each gene.
#ofNeighborsTxt_combined: Quadratic mean of number of neighbors in the literature-mined interaction network, associated with each gene.
#ofHighlyCoexpressed_combined: Quadratic mean of number of genes highly co-expressed with each gene.
#Common_PPI_Neighbors: Number of common neighbors in the PPI network for both genes.
#common_PWY_neighbors: Number of common neighbors in the pathways network for both genes.
#Common_Txt_Neighbors: Number of common neighbors in the literature-mined interaction network for both genes.
#Common_coexpressed: Number of common higly co-expressed genes for both genes.

The data files are available here.
  pathways
    kegg: gene to pathway codes python dictionary
    reactome: gene to pathway codes python dictionary
  phenotypes
    hpo: gene to phenotype codes python dictionary
  networks: PPI, pathways and literaure mined interaction network (dot) files and all shortest paths python dcitionary
  coex: mutual coexpression rank pythin dictionary
  evolgen: python dictionaries containing evolutionary and genomic features
    LOF intolerance
    haploinsufficiency
    dN/dS
    protein age
    gene essentiality
    
