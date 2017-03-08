library(WGCNA)

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
calculate_adjacency <- function(datan=NULL,filter=NULL,softPower=2) {
  return(sapply(1:nrow(filter),function(i){abs(cor(datan[,filter[i,1]],datan[,filter[i,2]]))^softPower}))
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



#' Run mcupgma
#' 
#' @param filter    Filter-Matrix
#' @param dist      Filter-Matrix
#' @param max_singleton     The maximal singleton in the clustering. Usually equals the number of features.
#' @return Dendrogram
#' @export
nb_mcupgma <- function(filter=NULL,dist=NULL,max_singleton=dim(datan)[2],cores=getOption("mc.cores", 2L)) {
  dir.create(path="clustering",recursive=TRUE)
  system(paste0("rm -rf clustering/*"))
  system(paste0("rm -r iteration_*/"))
  
  write.table(file="clustering/dist.edges",cbind(format(filter,scientific=FALSE,trim=TRUE),format(dist,scientific=FALSE,trim=TRUE)),row.names=FALSE,col.names=FALSE,sep="\t",quote=FALSE)
  system("gzip -f clustering/dist.edges")
  cluster <- file.path(netboost:::netboostPackagePath(), "mcupgma", "scripts", "cluster.pl")
  system(paste0(cluster," -max_distance ",1," -max_singleton ",max_singleton," -iterations 1000 -heap_size 10000000 -num_hash_buckets 40 -jobs ",cores," -retries 1 -output_tree_file clustering/dist.mcupgma_tree -split_unmodified_edges ",cores," clustering/dist.edges.gz")) 
  
  return(as.matrix(read.table(file="clustering/dist.mcupgma_tree",row.names=NULL,col.names=c("cluster_id1","cluster_id2","distance","cluster_id3"))))
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
  forest <- as.matrix(forest)
  
    if (is.null(forest) || !is.matrix(forest))
    stop("forest must be provided (as integer matrix)")
  
  return(netboost:::cpp_tree_search(forest))
}







#' Calculate dendrogram for an individual tree from a forest
#'
#' @param tree A list with two elements. ids, which is an integer vector of feature identifiers and rows, which is an integer vector of selected rows in the corresponding forest
#' @param datan     Dataset
#' @param forest Matrix
#' @return List of tree specific objects including dendrogram, tree data and features.
#' @export
tree_dendro <- function(tree=NULL, datan=NULL, forest=NULL) {
  forest <- as.data.frame(forest)
  index.features <- tree$ids[tree$ids <= dim(datan)[2]]
  data_tree <- datan[,index.features]
  
  colnames_tree <- colnames(datan)[index.features]
  tree_cluster <- forest[tree$rows,]

  #adjust for removed features
  all_ids <- rev(sort(unique(c(forest$cluster_id1,forest$cluster_id2,forest$cluster_id3))))
  none_tree_ids <- all_ids[all_ids %in% tree$ids]
  for(i in none_tree_ids){
    for(j in c(1,2,4)){
      tree_cluster[tree_cluster[,j]>i,j] <- tree_cluster[tree_cluster[,j]>i,j]-1
    }
  }

  feature_names <- colnames(data_tree)
  
  cutpoint <- dim(data_tree)[2]
  dendro <- list()
  dendro$merge <- tree_cluster[,c(1,2)]
  dendro$merge[dendro$merge<=cutpoint] <- -dendro$merge[dendro$merge<=cutpoint]
  dendro$merge[dendro$merge>0] <- dendro$merge[dendro$merge>0]-cutpoint
  dendro$merge <-  apply(dendro$merge,c(1,2),function(x){(as.integer(x))})
  dendro$height <- tree_cluster[,3]
  dendro$order <- 1:cutpoint
  dendro$labels <- colnames_tree
  class(dendro) <- "hclust"
  b <- as.dendrogram(dendro)
  b.order <- order.dendrogram(b)
  dendro$order <- b.order
 
  return(list(dendro=dendro,data=tree_data,names=feature_names)) 
}  
  
#' Module detection for a single tree
#'
#' @param tree_dendro List of tree specific objects including dendrogram, tree data and features originating from the tree_dendro function.
#' @param datan     Dataset
#' @param forest Matrix
#' @param MEDissThres Module Eigengene Dissimilarity Threshold for merging close modules.
#' @return Object of class hclust
#' @export
cut_dendro <- function(tree_dendro=NULL, minClusterSize= 10L, datan, MEDissThres = NULL) {
  
  dynamicMods = cutreeDynamic(dendro = tree_dendro$dendro, method="tree", deepSplit = TRUE, minClusterSize = minClusterSize)#, cutHeight=0.99);
  ### Merging of Dynamic Modules ###
  # Calculate eigengenes
  MEList = moduleEigengenes(tree_dendro$data, colors = dynamicMods)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  plot(METree, main = "Netboost: Clustering of module eigengenes",
       xlab = "", sub = "")
  abline(h=MEDissThres, col = "red")
  
  merged = mergeCloseModules(tree_dendro$data, dynamicMods, cutHeight = MEDissThres, verbose = 3)
  mergedColors = merged$colors;
  # Calculate eigengenes
  MEList = moduleEigengenes(tree_dendro$data, colors = merged$colors)
  MEs = MEList$eigengenes
  MEDiss = 1-cor(MEs);
  METree = hclust(as.dist(MEDiss), method = "average");
  plot(METree, main = "Netboost: Clustering of merged module eigengenes",
       xlab = "", sub = "")

  plotDendroAndColors(dendro=tree_dendro$dendro, mergedColors,"Merged Dynamic", dendroLabels = FALSE, hang = 0.01, addGuide = TRUE, guideHang = 0.05, main="Cluster Dendrogram Netboost")
  cat("Netboost extracted",length(table(mergedColors)[-1]),"modules with an average size of",mean(table(mergedColors)[-1])," from this tree.\n")
  
  return(mergedColors)
}

#' Cut trees
#'
#' @param trees List of trees, where one tree is a list of ids and rows
#' @return List
#' @export
cut_trees <- function(trees=NULL) {
  res <- list()
  i <- 1L
  for(tree in trees){
    pdf(file=paste0("tree_",i,".pdf"))
    res[[i]] <- cut_dendro(tree_dendro(tree))
    i <- i+1
    dev.off()
  }
  return(res)
}


#' TCGA RNA and methylation measurement on chromosome 18 for 180 AML patients.
#'
#' @format A data frame with 180 rows and 5283 variables:
#' @source \url{http://www.tcga.com/}
"tcga_aml_meth_rna_chr18"

#' TCGA covariates for 180 AML patients.
#'
#' @format A data frame with 188 rows and 78 variables:
#' @source \url{http://www.tcga.com/}
"tcga_aml_covariates"