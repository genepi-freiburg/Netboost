#' Test/example code.
#'
#' @export
exampleNB <- function() {
  # load data
  # methylation and RNA data
  data(tcga_aml_meth_rna_chr18)
  data(tcga_aml_covariates)

  filter <- nb_filter_boosting(datan=tcga_aml_meth_rna_chr18, stepno=20L, until=0L, progress=500L, cores=2L, mode=2L)

  adjacencies <- calculate_adjacency(datan=tcga_aml_meth_rna_chr18, filter=filter,softPower=6)

  dist <- dist_tom(filter=filter, adjacency=adjacencies, cores=2L)

  forest <- nb_mcupgma(filter=filter,dist=dist,max_singleton=dim(tcga_aml_meth_rna_chr18)[2],cores=2L)
 
  trees <- tree_search(forest);
  test <- tree_dendro(tree=trees[[1]], datan=tcga_aml_meth_rna_chr18, forest=forest)
  results <- cut_trees(trees)
  
  pdf(file="pdfs/dendrogram_netboost.pdf",height=10,width=500)
  plot(dendro.netboost,xlab='',sub='',lab=FALSE,hang=0.01)
  dev.off()
  
  
}