# Package internal variables (as environment).
.netboostInternal <- new.env(parent = emptyenv())

#' Package startup: used to fetch installation path of the own package,
#' as required for executing binary programs delivered with it.
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
