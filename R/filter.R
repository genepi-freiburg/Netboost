library(parallel)

#' Boosting via C++ function. Parallelisation by R-package parallel with forking
#' (overhead of this method does not fall into account as single steps are ~10s).
#'
#' Parallelisation via multicore (via 'parallel'-package). So *nix only atm.
#'
#' @param datan     Dataset
#' @param stepno    Integer amount of boosting steps
#' @param until     Stop at index/column (if 0: iterate through all columns)
#' @param progress  Integer. If > 0, print progress after X steps (mind: parallel!)
#' @param cores     Integer. Amount of CPU cores used (<=1 : sequential)
#' @param mode      Integer. Mode (0: x86, 1: FMA, 2: AVX). Features are only available
#'                  if compiled accordingly and available on the hardware.
#' @return List (for each column vector with items received by boosting)
#' @export
nb_filter_boosting <- function(datan, stepno=20, until=0,
                               progress=100,
                               cores=getOption("mc.cores", 2L),
                               mode=2) {
  if (until == 0)
    until = ncol(datan);
  
  ## Initialize data structures for optimized boosting (once)
  netboost:::cpp_filter_base(datan, stepno, mode=mode);
  
  ## Parallelization "conventional" via mclapply. Not really accountable overhead,
  ## as single calls take ~10 seconds.
  if (cores > 1) {
    print(paste("Parallel version:", cores, "cores"))
    
    ret <- mclapply(seq(1, until),
                    function(x) {
                      if ((progress > 0) && (((x-1) %% progress) == 0)) {
                        print(sprintf("idx: %d (%.1f%%) - %s", x, x * 100 / until, date()))
                      }
                      
                      netboost:::cpp_filter_step(x)
                    },
                    mc.cores=cores)
  } else {   ## Sequential function for debugging.
    print(paste("Sequential version"))
    
    ret <- lapply(seq(1, until),
                  function(x) {
                    if (progress && (((x-1) %% progress) == 0)) {
                      print(timestamp(quiet=TRUE))
                      print(sprintf("idx: %d (%.1f%%)", x, x * 100 / until))
                    }
                    
                    netboost:::cpp_filter_step(x)
                  })
  }
  
  ## Important!: stop (free memory, else suitable memory is still blocked)
  netboost:::cpp_filter_end();
  
  return(ret)
}
