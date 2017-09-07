#' Test/example code.
#'
#' @export
nb_example <- function() {
  # load data
  # methylation and RNA data
  data(tcga_aml_meth_rna_chr18) # 180 patients x 5283 features
  data(tcga_aml_covariates)

  options("mc.cores"=10L)
  
  filter <- nb_filter(datan=tcga_aml_meth_rna_chr18, stepno=20L)

  dist <- nb_dist(datan=tcga_aml_meth_rna_chr18, filter=filter, softPower=6)

  pdf(file="results_netboost.pdf",width = 30)
    results <- nb_clust(datan=tcga_aml_meth_rna_chr18, filter=filter, dist=dist, minClusterSize = 10L, MEDissThres = 0.25)
  dev.off()
  plot_trees(results)
  
 
  
  ### compute cluster settings
  # install.packages("/data/bin/netboost_1.021.tar.gz", repos = NULL, type="source")
  
  # library(netboost)

  
  system(paste0("rm -rf clustering/*"))
  system(paste0("rm -r clustering/"))
  system(paste0("rm -r iteration_*/"))
  
  
  
}