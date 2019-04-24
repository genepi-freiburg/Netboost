## Note: here also the setup NAMESPACE directives are set (as Roxygen is
## creating the namespace, those are given as Roxygen attributes). The
## importFrom statements are required to instantly load the linked libraries
## Rcpp and RcppParallel, else loading of the own shared lib would fail (the
## imported functions do not matter, but for each package at least one import
## must be present).

## Those are required for CRAN checks but bug out on BiocCheck
#' Package startup: used to fetch installation path of the own package,
#' as required for executing binary programs delivered with it.
#' 
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel setThreadOptions
#' @importFrom R.utils gzip isGzipped
#' @importFrom parallel mclapply
#' @importFrom colorspace rainbow_hcl
#' @importFrom grDevices dev.off gray pdf
#' @importFrom graphics abline layout par plot
#' @importFrom stats as.dendrogram as.dist cor cov prcomp hclust order.dendrogram
#' @importFrom dynamicTreeCut cutreeDynamic indentSpaces printFlush
#' @importFrom impute impute.knn
#' @importFrom WGCNA allowWGCNAThreads mergeCloseModules plotDendroAndColors
#' @importFrom WGCNA moduleColor.getMEprefix pickSoftThreshold
#' @importFrom utils data packageDescription read.table write.table
#' @importFrom methods is
#'
#' @useDynLib netboost
#'
#' @examples 
#' \dontrun{nb_example()}
#' @return none
#' @param libname Path to R installation (base package dir)
#' @param pkgname Package name (should be "netboost")
.onAttach <- function(libname, pkgname) {
    desc <- packageDescription(pkgname)

    # If no default core count given, detect.  
    if (is.null(getOption("mc.cores")) || !is.integer(getOption("mc.cores"))) {
        # logical = FALSE is not working correctly if CPU has logical cores, which
        # are disabled (at least Linux).
        # Means: if CPU has logical cores, core count should be set manually.
        # cores <- parallel::detectCores() # Bioconductor does not like this
        cores <- NA

        if (is.na(cores)) cores <- 1

        options("mc.cores" = cores)
    }
    
    ## Optional startup message, mainly for development.
    packageStartupMessage(
        paste(pkgname,
              desc$Version,
              "loaded"),
        paste(
            "Default CPU cores:",
            getOption("mc.cores"),
            "\n"),
        appendLF = TRUE
    )
    #                              "Loaded from:", libname),

    # Create temp subfolder in tempdir()
    netboostTmpCleanup()

    ## Add the current (real) loading path to MCUPGMA Makefiles
    ## (install_path.mk is loaded by definitions.mk, which is
    ## included in all real Makefiles).
    mcupgma_install <- file.path(netboostMCUPGMAPath(),
                                 "install_path.mk")

    ## If this file is not existing in this location, this is a non working
    ## installation (may happen during build and included test-loads) (writeLines
    ## throws warning in R CMD check, but we do valid stuff here)
    if (file.exists(mcupgma_install)) {
        # R complains about writeLines (false positive, as not writing to STDOUT).
        # Replaced with write.table to pass package check.
        txt <- c(paste("export INSTALL_PATH := ", netboostMCUPGMAPath()),
                 paste("export TMP_PATH := ", netboostTmpPath()))
        write.table(file = mcupgma_install, as.data.frame(txt),
                    quote = FALSE, row.names = FALSE,
                    col.names = FALSE, append = FALSE, sep="")
        ##    filew <-file(mcupgma_install, open="w")
        ##    writeLines(con=filew, text=c(paste("export INSTALL_PATH := ",
        ##    mcupgmaPath)))
        ##    writeLines(con=filew, text=c(paste("export TMP_PATH := ",
        ##    netboostTmpPath())))
        ##    close(filew)
    }
    # Else successful build would be warned.
    else {
        warning(paste("File not written:", mcupgma_install,
                      "(okay during build, error after installation)"))
    }
}

## #' If package detached, clean up temporary folders.
## #' @return none
## #' @param libpath Library path (unused)
##.onDetach <- function(libpath) {
##    print("kthnxbye")
##}

#' Returns the absolute path to "exec" folder in the package.
#'
#' @return Absolute path of installed package
netboostPackagePath <- function() {
    return(system.file(package="netboost"))
}

#' Returns the absolute path to temporary folder of the package.
#' To change temporary path, use normal R variables (TEMPDIR etc).
#'
#' @return Absolute path for "exec" folder
netboostTmpPath <- function() {
    return(file.path(tempdir(), "netboost"))
}

#' Returns the absolute path to folder with mcupgma executables and scripts.
#'
#' @return Absolute path for "mcupgma" folder
netboostMCUPGMAPath <- function() {
    return(file.path(netboostPackagePath(), "mcupgma"))
}
