#' Execute a program/script from the installes MCUPGMA suite.
#'
#' @param exec    Name of the file of the executable.
#' @param console Print output to R console or fetch for return to caller
#' @return console=TRUE: exit code (0: no error). console=FALSE: STDOUT/STDERR output
#' @export
mcupgma_exec <- function(exec=NULL, ..., console=TRUE) {
  if (is.null(exec)) stop("mcupgma_exec: call without executable")
  
  exec_abs_path <- file.path(netboostMCUPGMAPath(), exec)
  
  if (!file.exists(exec_abs_path))
    stop(paste("mcupgma_exec: file does not exist:", exec, exec_abs_path))

  # Types may become mixed giving many parameters, so ensure boolean.
  console <- isTRUE(as.logical(console))

  # Arguments must be bundled in character vector
  call_args <- as.character(unlist(list(...), use.names=FALSE))
  
  # Default: STDOUT/STDERR to R console
  std <- ""
  
  # If console is not set, write output to temporary file.
  if (!console) std <- tempfile()
  
  # Execute command (always fetch return code)
#  print(paste("Starting:", exec_abs_path))
#  print(paste(..., sep=",", collate=""))
#  flush.console()
  
  ret <- system2(command = exec_abs_path,
                 stdout = std, stderr = std,
                 wait=TRUE,
                 args=call_args,
                 env=c(paste0("TMP_PATH=", netboostTmpPath())))

#  flush.console()
  
  # If execution returned error, always throw warning.
  if (ret > 0)
    warning(paste("Execution of", exec, "returned error:", ret),
            call.=FALSE)
  
  # If STDIN/STDERR was printed on console, return exit code.
  # Else, read output from temporary file and return as text.
  if (console) {
    return(ret)
  } else {
    std_output <- readChar(std, file.info(std)$size)
    unlink(std)
    return(std_output)
  }
}
