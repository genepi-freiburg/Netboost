# Package internal variables (as environment).
.netboostInternal <- new.env(parent = emptyenv())

# Note: here also the setup NAMESPACE directives are set (as Roxygen is creating
# the namespace, those are given as Roxygen attributes).
# The importFrom statements are required to instantly load the linked libraries
# Rcpp and RcppParallel, else loading of the own shared lib would fail (the imported
# functions do not matter, but for each package at least one import must be present).

#' Package startup: used to fetch installation path of the own package,
#' as required for executing binary programs delivered with it.
#' 
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel setThreadOptions
#' @useDynLib netboost
#'
#' @param libname Path to R installation (base package dir)
#' @param pkgname Package name (should be "netboost")
.onAttach <- function(libname, pkgname) {
  desc <- packageDescription(pkgname)
  
  ## Optional startup message, mainly for development.
  packageStartupMessage(paste(pkgname,
                              desc$Version,
                              desc$Date,
                              "Loaded from:", libname),
                        appendLF = TRUE)

  ## Path to "exec"-folder in installed package
  pPath <- file.path(libname, pkgname)

  ## Store exec and general path in package variables
  assign("exec_path", file.path(pPath, "exec"), envir = .netboostInternal)
  assign("mcupgma_path", file.path(pPath, "mcupgma"), envir = .netboostInternal)
  assign("pkg_path", pPath, envir = .netboostInternal)
}

#' Returns the absolute path to "exec" folder in the package.
#'
#' @return Absolute path for "exec" folder
netboostExecPath <- function() {
  if (exists("exec_path", envir = .netboostInternal)) {
    return(get("exec_path", envir = .netboostInternal))
  }
  else {
    stop("Executable path not existing (key exec_path missing in envir)")
  }
}

#' Returns the absolute path to folder with mcupgma executables and scripts.
#'
#' @return Absolute path for "mcupgma" folder
netboostMCUPGMAPath <- function() {
  if (exists("mcupgma_path", envir = .netboostInternal)) {
    return(get("mcupgma_path", envir = .netboostInternal))
  }
  else {
    stop("mcupgma path not existing (key mcupgma_path missing in envir)")
  }
}

#' Returns the absolute path to "exec" folder in the package.
#'
#' @return Absolute path of installed package
netboostPackagePath <- function() {
  if (exists("pkg_path", envir = .netboostInternal)) {
    return(get("pkg_path", envir = .netboostInternal))
  }
  else {
    stop("Package path not existing (key pkg_path missing in envir)")
  }
}
