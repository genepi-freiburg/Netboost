#' Example to get access to the MCUPGMA executables.
#'
#' @export
mcupgma_example <- function() {
  exec <- netboostMCUPGMAPath()
  files <- Sys.glob(file.path(exec, '*'))
  paste("Available MCUPGMA executables and scripts under:", exec)
  print(sapply(files, basename, USE.NAMES=FALSE))
}


#' Test/example code.
#'
#' @param cores CPU cores to use
#' @export
nb_example <- function(cores=2L) {
  ## Supress warning of third party package (which one?)
  if (cores > 1) allowWGCNAThreads()

  # load data
  # methylation and RNA data
  data(tcga_aml_meth_rna_chr18) # 180 patients x 5283 features
  data(tcga_aml_covariates)
  options("mc.cores"=cores)

  pdf(file=file.path(getwd(), "results_netboost.pdf"), width = 30)
  results <- netboost(datan=tcga_aml_meth_rna_chr18,stepno=20L, minClusterSize = 10L, MEDissThres = 0.25) 
  dev.off()
  
  ### Transfer results to the same data (bug check)
  tmp <-  nb_transfer(nb_summary = results, new_data = tcga_aml_meth_rna_chr18)
  sum(results$MEs!=tmp)
  
  ### compute cluster settings
  # install.packages("/data/bin/netboost_1.023.tar.gz", repos = NULL, type="source")
  # install.packages("~/git/netboost_1.023.tar.gz", repos = NULL, type="source")
  
  # library(netboost)

  # Debugging: check produced temporary files.
  stop("Check files")

  #  system(paste0("rm -rf clustering/*"))
  #  system(paste0("rm -r clustering/"))
  #  system(paste0("rm -r iteration_*/"))
  
  # Cleanup all produced temporary filed (esp. clustering/iteration_*)
  netboostTmpCleanup()
}
