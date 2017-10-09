#' Example to get access to the MCUPGMA executables.
#'
#' @export
mcupgma_example <- function() {
  exec <- netboostMCUPGMAPath()
  files <- Sys.glob(file.path(exec, '*'))
  paste("Available MCUPGMA executables and scripts under:", exec)
  tmp <- sapply(files, basename)
  names(tmp) <- NULL
  print(tmp)
}

#' Execute a program/script from the installes MCUPGMA suite.
#'
#' @param exec  Name of the file of the executable.
#' @return response from system2() call
#' @export
mcupgma_exec <- function(exec=NULL, ...) {
  if (is.null(exec)) stop("mcupgma_exec: call without executable")

  exec_abs_path <- file.path(netboostMCUPGMAPath(), exec)

  if (!file.exists(exec_abs_path))
    stop(paste("mcupgma_exec: file does not exist:", exec))

  return(system2(command=exec_abs_path,
                 args=...,
                 stdout=TRUE, stderr=TRUE))
}

#' Test/example code.
#'
#' @export
nb_example <- function() {
  # load data
  # methylation and RNA data
  data(tcga_aml_meth_rna_chr18) # 180 patients x 5283 features
  data(tcga_aml_covariates)

  options("mc.cores"=20L)
  options("mc.cores"=8L)
  
  filter <- nb_filter(datan=tcga_aml_meth_rna_chr18, stepno=20L)
  dist <- nb_dist(datan=tcga_aml_meth_rna_chr18, filter=filter, softPower=6)
  pdf(file="results_netboost.pdf",width = 30)
    results <- nb_clust(datan=tcga_aml_meth_rna_chr18, filter=filter, dist=dist, minClusterSize = 10L, MEDissThres = 0.25)
  dev.off()
  
  pdf(file="results_netboost.pdf",width = 30)
  results <- netboost(datan=tcga_aml_meth_rna_chr18,stepno=20L, softPower=6L, minClusterSize = 10L, MEDissThres = 0.25) 
  dev.off()
  
  ### Transfer results to the same data (bug check)
  tmp <-  nb_transfer(nb_summary = results, new_data = tcga_aml_meth_rna_chr18)
  sum(results$MEs!=tmp)
  
  ### compute cluster settings
  # install.packages("/data/bin/netboost_1.023.tar.gz", repos = NULL, type="source")
  # install.packages("~/git/netboost_1.023.tar.gz", repos = NULL, type="source")
  
  # library(netboost)

  system(paste0("rm -rf clustering/*"))
  system(paste0("rm -r clustering/"))
  system(paste0("rm -r iteration_*/"))
}