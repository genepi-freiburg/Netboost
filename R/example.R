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
#' @param cores Integer. CPU cores to use.
#' @param keep Logical. Keep mcupgma intermediate files.
#' @export
nb_example <- function(cores = getOption("mc.cores", 2L),
                       keep = FALSE) {
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
  
  # Cleanup all produced temporary filed (esp. clustering/iteration_*)
  if (!keep) netboostTmpCleanup()
}
