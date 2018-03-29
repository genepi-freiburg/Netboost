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
#' @importFrom grDevices dev.off gray pdf rainbow
#' @useDynLib netboost
#'
#' @param libname Path to R installation (base package dir)
#' @param pkgname Package name (should be "netboost")
.onAttach <- function(libname, pkgname) {
  desc <- packageDescription(pkgname)
  
  ## Optional startup message, mainly for development.
  packageStartupMessage(paste(pkgname,
                              desc$Version,
                              "loaded"),
#                              desc$Date,
#                              "Loaded from:", libname),
                        appendLF = TRUE)

  ## Path to "exec"-folder in installed package
  pPath <- file.path(libname, pkgname)

  mcupgmaPath <- file.path(pPath, "mcupgma")
  
  ## Store exec and general path in package variables
  assign("exec_path", file.path(pPath, "exec"), envir = .netboostInternal)
  assign("mcupgma_path", mcupgmaPath, envir = .netboostInternal)
  assign("pkg_path", pPath, envir = .netboostInternal)

  ## Per default, temporary data is written to R tempdir. But as those may become
  ## large and R normally use /tmp, user must be enabled to change those later.
  nb_set_tempdir(file.path(tempdir(), "netboost"))

  ## Add the current (real) loading path to MCUPGMA Makefiles (install_path.mk
  ## is loaded by definitions.mk, which is included in all real Makefiles).
  mcupgma_install <- file.path(mcupgmaPath, "install_path.mk")

  ## If this file is not existing in this location, this is no working installation
  ## (may happen during build and included test-loads)
  if (file.exists(mcupgma_install)) {
#    print(paste("Modified:", mcupgma_install))
    filew <-file(mcupgma_install, open="w")
    writeLines(con=filew, text=c(paste("export INSTALL_PATH := ", mcupgmaPath)))
    writeLines(con=filew, text=c(paste("export TMP_PATH := ", netboostTmpPath())))
    close(filew)
  }
  else {
    warning(paste("File not written (as it does not exist):", mcupgma_install))
  }
}

#' If package detached, clean up temporary folders.
.onDetach <- function(libpath) {
  netboostTmpCleanup()
}

#' Assigns temporary path for internal use (esp. mcupgma)
#' 
#' @param tmp Directory (Default: R temporary folder)
#' @export
nb_set_tempdir <- function(tmp = NULL) {
  # TODO Cleanup maybe currently existing temporary folder.
  netboostTmpCleanup()

  folder <- tmp

  if (is.null(tmp))
    folder <- tempdir()
  
  if (file.exists(folder) && !dir.exists(folder))
    stop(paste("Given temporary exists as file:", folder),
         call.=FALSE)
  
  if (!dir.exists(folder)) {
    if (is.null(tmp))
      message("Using temporary directory:", folder)
    else
      warning(paste("Given temporary directory not existing. Created:", folder),
              call. = FALSE)
    
    if (!dir.create(folder, recursive = TRUE, showWarnings = TRUE))
      stop(paste("Error creating temporary folder:", folder))
  }
  else {
    # Existing folder must be empty, as will be removed on cleanup.
    files <- Sys.glob(file.path(folder, '*'))
    
    if (length(files) > 0)
      stop(paste("Given temporary folder not empty:", folder))
  }

  assign("temp_dir", folder, envir = .netboostInternal)
}

#' Cleans the netboost temporary folder.
netboostTmpCleanup <- function() {
  folder <- netboostTmpPath(nostop=TRUE)
  
  if (folder == "") return()

  if (dir.exists(folder)) {
    message(paste("Cleaning temporary folder:", folder))
    
    ## Delete and recreate more convenient than globbing through the folders
    unlink(folder, recursive = TRUE)
    dir.create(folder)
  }
}

#' Returns the absolute path to "exec" folder in the package.
#'
#' @param  nostop   Return on error (default: stop)
#' @return Absolute path for "exec" folder
netboostTmpPath <- function(nostop=FALSE) {
  if (exists("temp_dir", envir = .netboostInternal)) {
    return(get("temp_dir", envir = .netboostInternal))
  }
  else {
    if (nostop)
      return("")
    else
      stop("No temporary folder set")
  }
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
