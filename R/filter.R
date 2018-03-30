library(parallel)

#' Boosting via C++ function. Parallelisation by R-package parallel with forking
#' (overhead of this method does not fall into account as single steps are ~10s).
#'
#' Parallelisation via multicore (via 'parallel'-package). So *nix only atm.
#'
#' @param datan     Data frame were rows correspond to samples and columns to features.
#' @param stepno    Integer amount of boosting steps
#' @param until     Stop at index/column (if 0: iterate through all columns)
#' @param progress  Integer. If > 0, print progress after every X steps (mind: parallel!)
#' @param cores     Integer. Amount of CPU cores used (<=1 : sequential)
#' @param mode      Integer. Mode (0: x86, 1: FMA, 2: AVX). Features are only available
#'                  if compiled accordingly and available on the hardware.
#' @return matrix n times 2 matrix with the indicies of the n unique entrees of the filter
#' @export
nb_filter <- function(datan, stepno=20L, until=0L,
                      progress=1000L,
                      cores=getOption("mc.cores", 2L),
                      mode=2L) {
  if (!exists("datan"))
    stop("datan must be provided")
  
  if (!(is.data.frame(datan) && (nrow(datan) > 0) && (ncol(datan) > 0)))
    stop("datan must be a data frame with dim() > (0,0).")
  
  if (!(is.integer(stepno) && (stepno > 0)))
    stop("stepno must be an integer > 0.")
  
  if (!(is.integer(until) && (until >= 0)))
    stop("until must be an integer >= 0.")
  
  if (until == 0)
    until = ncol(datan);
  
  # check mode
  if(!(mode %in% c(0, 1, 2))){
    stop("mode must be 0 (x86), 1 (FMA) or 2 (AVX).")
  }

  ## Initialize data structures for optimized boosting (once)
  cpp_filter_base(as.matrix(datan), stepno, mode=mode);
  
  ## Parallelization "conventional" via mclapply.
  if (cores > 1) {
    #print(paste("Parallel version:", cores, "cores"))
    
    boosting_filter <- mclapply(seq(1, until),
                    function(x) {
                      if ((((x-1) %% progress) == 0)) {
                        print(sprintf("idx: %d (%.1f%%) - %s", x, x * 100 / until, date()))
                      }
                      
                      cpp_filter_step(x)
                    },
                    mc.cores=cores)
  } else {   ## Sequential function for debugging.
#    print(paste("Sequential version"))
    boosting_filter <- lapply(seq(1, until),
                              function(x) {
                                if ((((x-1) %% progress) == 0)) {
                                  print(sprintf("idx: %d (%.1f%%) - %s",
                                                x,
                                                x * 100 / until,
                                                date()))
                                }
                                
                                cpp_filter_step(x)
                              })
  }
  
  ## Important!: stop (free memory, else suitable memory is still blocked)
  cpp_filter_end();

  filter <- do.call("rbind", lapply(1:length(boosting_filter), function(x) {
    return(as.data.frame(cbind(as.integer(boosting_filter[[x]]),
                               as.integer(rep(x,length(boosting_filter[[x]]))))))
  }))
  
  filter <- unique(t(apply(filter,1,sort)))
  colnames(filter) <- c("cluster_id1","cluster_id2")
  rownames(filter) <- 1:nrow(filter)
  
  return(filter)
}