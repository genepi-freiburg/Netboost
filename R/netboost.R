#' Skeleton (Test)
#'
#' @return Nix
#' @export
netboost <- function() {
  print("Netboost. Yeah")
  
  prg <- file.path(netboostExecPath(), "test.pl")

  system(prg)
  
  prg <- file.path(netboostPackagePath(), "mcupgma", "scripts", "cluster.pl")
  
  system(prg)
  
  
}


#' Calculate network adjacencies for filter
#' 
#' @param datan     Dataset
#' @param filter    Filter-Matrix
#' @param softPower Integer. Exponent of the transformation
#' @return Vector with adjacencies for the filter
#' @export
calculate_adjacency <- function(filter=NULL,softPower=1) {
  adjacency <- sapply(1:nrow(filter),function(i){abs(cor(datan[,filter[i,1]],datan[,filter[i,2]]))^softPower})
  return(adjacency)
}

#' Calculate distance
#' (external wrapper for internal C++ function)
#' Parallelisation inside C++ program with RcppParallel.
#' 
#' @param filter Filter-Matrix
#' @param adjacency Adjacency-Vector
#' @param Integer. Amount of CPU cores used (<=1 : sequential)
#' @return Vector with distances (same length as adjacency)
#' @export
dist_tom <- function(filter=NULL,
                     adjacency=NULL,
                     cores=getOption("mc.cores", 2L)) {
  if (is.null(filter) || is.null(adjacency))
    stop("Both filter and adjacency must be provided")

  if (!(is.matrix(filter) && (nrow(filter) > 0) && (ncol(filter) > 0)))
    stop("filter must be matrix with dim() > (0,0)")

  if (!(is.vector(adjacency) && (length(adjacency) > 0)))
    stop("adjacency is required a vector with length > 0")

  cores <- max(cores, 1)

  ## RcppParallel amount of threads started  
  setThreadOptions(numThreads=cores)

  return(netboost:::cpp_dist_tom(filter, adjacency))
}

#' Tree search
#' (external wrapper for internal C++ function)
#'
#' @param forest Matrix
#' @return List
#' @export
tree_search <- function(forest=NULL) {
  # Check for integer values cannot be done here, as either the user must
  # have set up all values with as.integer() or R delivers default numeric
  # (double). In that case, Rcpp converts the matrix.
  # (Alternative: convert manually "matrix(as.integer(forest), nrow=nrow(forest))")
  if (is.null(forest) || !is.matrix(forest))
    stop("forest must be provided (as integer matrix)")

  return(netboost:::cpp_tree_search(forest))
}
