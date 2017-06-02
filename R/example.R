#' Test/example code.
#'
#' @export
exampleNB <- function() {
  # load data
  # methylation and RNA data
  data(tcga_aml_meth_rna_chr18) # 180 patients x 5283 features
  data(tcga_aml_covariates)

  filter <- nb_filter_boosting(datan=tcga_aml_meth_rna_chr18, stepno=20L, until=0L, progress=500L, cores=2L, mode=2L)

  adjacencies <- calculate_adjacency(datan=tcga_aml_meth_rna_chr18, filter=filter,softPower=6)

  dist <- dist_tom(filter=filter, adjacency=adjacencies, cores=2L)

  forest <- nb_mcupgma(filter=filter,dist=dist,max_singleton=dim(tcga_aml_meth_rna_chr18)[2],cores=2L)
 
  trees <- tree_search(forest);

  results <- cut_trees(trees=trees,datan=tcga_aml_meth_rna_chr18, forest=forest, minClusterSize = 10L, MEDissThres = 0.25)

  plot_trees(results)
  
 
  
  ### compute cluster settings
  # library(netboost)
  data(tcga_aml_meth_rna_chr18)
  data(tcga_aml_covariates)
  filter <- nb_filter_boosting(datan=tcga_aml_meth_rna_chr18, stepno=20L, until=0L, progress=500L, cores=12L, mode=2L)
  adjacencies <- calculate_adjacency(datan=tcga_aml_meth_rna_chr18, filter=filter,softPower=6)
  dist <- dist_tom(filter=filter, adjacency=adjacencies, cores=12L)
  forest <- nb_mcupgma(filter=filter,dist=dist,max_singleton=dim(tcga_aml_meth_rna_chr18)[2],cores=12L)
  trees <- tree_search(forest)
  pdf(file="results_individual_trees.pdf")
  results <- cut_trees(trees=trees,datan=tcga_aml_meth_rna_chr18, forest=forest, minClusterSize = 10L, MEDissThres = 0.25)
  dev.off()
  pdf(file="results_netboost.pdf",width = 30)
  sum_res <- nb_summary(results, MEDissThres = 0.25)
  dev.off()
  
  
  system(paste0("rm -rf clustering/*"))
  system(paste0("rm -r clustering/"))
  system(paste0("rm -r iteration_*/"))
  
  
  
}