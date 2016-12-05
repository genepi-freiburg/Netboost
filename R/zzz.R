# Package internal variables (as environment).
.netboostInternal <- new.env()

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
  prg <- file.path(libname, pkgname, "exec")

  ## Store path in package variables
  assign("exec", prg, envir = .netboostInternal)
}

#' Call external program delivered in the package
#'
#' @return STDOUT from called program
#' @export
call_external <- function() {
  if (exists("exec", envir = .netboostInternal)) {
    prg <- file.path(get("exec", envir = .netboostInternal), "test.pl")
    return(system(command = prg, inter = TRUE))
  }
  else {
    stop("exec not given (key exec missing in envir)")
  }
}
