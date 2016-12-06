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
