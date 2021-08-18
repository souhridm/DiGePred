The features used to characterize each gene pair to train DiGePred were:

LoFintolerance: Quadratic mean of genes' sensitivity to loss-of-function variants. (https://www.deciphergenomics.org)
Haploinsufficiency: Quadratic mean of genes' resistance to haploinsufficiency. (https://www.deciphergenomics.org)
Age: Quadratic mean of gene/protein ages. (https://proteinhistorian.docpollard.org)
dN/dS: Quadratic mean of selection pressure on genes. (http://www.h-invitational.jp/evola/search.html)
Essentiality: Quadratic mean of organism's tolerance to loss of genes' funtionality. (https://v3.ogee.info/#/home)

The values of these features for individual genes are avaialble from the python dictionaries.
