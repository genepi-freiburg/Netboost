#' Example to get access to the MCUPGMA executables.
#'
#' @export
mcupgma_example <- function() {
  exec <- netboostMCUPGMAPath()
  files <- Sys.glob(file.path(exec, '*'))
  paste("Available MCUPGMA executables and scripts under:", exec)
  print(sapply(files, basename, USE.NAMES=FALSE))
}


#' Test/example code. Applies netboost to the TCGA-AML CHR18 DNA methylation and gene expression data supplied with the package. 
#'
#' @param cores Integer. CPU cores to use.
#' @param keep Logical. Keep mcupgma intermediate files.
#' @export
nb_example <- function(cores = getOption("mc.cores", 2L),
                       keep = FALSE) {
  # load data
  # methylation and RNA data
  data("tcga_aml_meth_rna_chr18",  # 180 patients x 5283 features
       package="netboost",
       envir = environment())

  pdfFile = file.path(tempdir(), "results_netboost.pdf")
  
#  pdf(file=file.path(getwd(), "results_netboost.pdf"), width = 30)
  pdf(file=pdfFile, width = 30)
  results <- netboost(datan = tcga_aml_meth_rna_chr18,
                      stepno = 20L,
                      softPower = 3L,
                      minClusterSize = 10L,
                      MEDissThres = 0.25)
  nb_plot_dendro(nb_summary = results,labels=TRUE,colorsrandom=TRUE,seed=123)
  dev.off()
  
  if (file.exists(pdfFile)) {
    message(paste0("PDF created:", pdfFile))
    
    # If default PDF viewer is assigned, try to show PDF.
    if (!is.null(getOption("pdfviewer"))) {
      system2(getOption("pdfviewer"), pdfFile)
    }
  }
  
  ### Transfer results to the same data (bug check)
  tmp <-  nb_transfer(nb_summary = results, new_data = tcga_aml_meth_rna_chr18)
  
  ## Check transfer of data
  sum(round(results$MEs,12) != round(tmp,12))
  sum(round(results$MEs,12) == round(tmp,12))

  # Cleanup all produced temporary filed (esp. clustering/iteration_*)
  if (!keep)
    netboostTmpCleanup()
  else
    print(paste("Kept MCUPGMA temporary files in:", netboostTmpPath()))
  
}


