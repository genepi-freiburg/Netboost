#' Skeleton.
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

#' Calculate distance
#' (external wrapper for internal C++ function)
#' 
#' @param filter Filter-Matrix
#' @param adjacency Adjacency-Vector
#' @return Vector with distances (same length as adjacency)
#' @export
dist_tom <- function(filter=NULL, adjacency=NULL) {
  if (is.null(filter) || is.null(adjacency))
    stop("Both filter and adjacency must be provided")

  if (!(is.matrix(filter) && (nrow(filter) > 0) && (ncol(filter) > 0)))
    stop("filter must be matrix with dim() > (0,0)")

  if (!(is.vector(adjacency) && (length(adjacency) > 0)))
    stop("adjacency is required a vector with length > 0")

  return(netboost:::rcpp_dist_tom(filter, adjacency))
}
