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
nb_mcupgma <- function(filter=NULL,dist=NULL,max_singleton=NULL,cores=getOption("mc.cores", 2L)) {
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
#' @param datan Dataset
#' @param forest Matrix
#' @return List of tree specific objects including dendrogram, tree data and features.
#' @export
tree_dendro <- function(tree=NULL, datan=NULL, forest=NULL) {
  index.features <- tree$ids[tree$ids <= dim(datan)[2]]
  data_tree <- datan[,index.features]
  
  colnames_tree <- colnames(datan)[index.features]
  tree_cluster <- forest[tree$rows,]

  all_ids <- rev(sort(unique(c(forest[,c(1,2,4)]))))
  none_tree_ids <- all_ids[!(all_ids %in% tree$ids)]
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
 
  return(list(dendro=dendro,data=data_tree,names=feature_names)) 
}  
  
#' Module detection for a single tree
#'
#' @param tree_dendro List of tree specific objects including dendrogram, tree data and features originating from the tree_dendro function.
#' @param datan Dataset
#' @param forest Matrix
#' @param MEDissThres Module Eigengene Dissimilarity Threshold for merging close modules.
#' @return Object of class hclust
#' @export
cut_dendro <- function(tree_dendro=NULL, minClusterSize= 10L, datan=NULL, MEDissThres = NULL, name_of_tree="", plot = TRUE) {
  dynamicMods <- cutreeDynamic(dendro = tree_dendro$dendro, method="tree", deepSplit = TRUE, minClusterSize = minClusterSize)
  ### Merging of Dynamic Modules ###
  # Calculate eigengenes
  MEList <- moduleEigengenes(tree_dendro$data, colors = dynamicMods)
  MEs <- MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss <- 1-cor(MEs);
  # Cluster module eigengenes
  if(length(MEDiss) > 1){
    METree <- hclust(as.dist(MEDiss), method = "average");
    if(plot == TRUE){
      plot(METree, main = paste0(name_of_tree,"Clustering of module eigengenes"),
           xlab = "", sub = "")
      abline(h=MEDissThres, col = "red")
    }
    
    merged <- mergeCloseModules(exprData=tree_dendro$data, dynamicMods, cutHeight = MEDissThres, verbose = 3)
    mergedColors <- merged$colors;
    # Calculate eigengenes
    MEList <- moduleEigengenes(tree_dendro$data, colors = merged$colors)
    MEs <- MEList$eigengenes
    MEDiss <- 1-cor(MEs);
    METree <- hclust(as.dist(MEDiss), method = "average");
    if(plot == TRUE){
      plot(METree, main =paste0(name_of_tree,"Clustering of merged module eigengenes"),
         xlab = "", sub = "")
      plotDendroAndColors(dendro=tree_dendro$dendro, colors=mergedColors,"Merged Dynamic", dendroLabels = FALSE, hang = 0.01, addGuide = TRUE, guideHang = 0.05, main=paste0(name_of_tree,"Cluster Dendrogram"))
    }
    }else{
    cat("\nOnly one module in ",name_of_tree,".\n")
    mergedColors <- dynamicMods
    if(plot == TRUE){
      plot(tree_dendro$dendro, main=paste0(name_of_tree,"Cluster Dendrogram (Tree maximaly consists out of one module.)"))
    }
  }
  cat("\nNetboost extracted",length(table(mergedColors)),"modules (including background) with an average size of",mean(table(mergedColors)[-1])," (excluding background) from ",substr(name_of_tree,start=1,stop=(nchar(name_of_tree)-1)),".\n")
  return(list(colors=mergedColors,MEs=MEs))
}

#' Cut trees
#'
#' @param trees List of trees, where one tree is a list of ids and rows
#' @return List
#' @export
cut_trees <- function(trees=NULL, datan=NULL, forest=NULL, minClusterSize= 10L, MEDissThres = NULL, plot = TRUE) {
  res <- list()
  i <- 1L
  for(tree in trees){
    tree_dendro <- tree_dendro(tree=tree,datan=datan,forest=forest)
    res[[i]] <- list()
    res[[i]][["dendro"]] <- tree_dendro$dendro
    res[[i]][["data"]] <- tree_dendro$data
    res[[i]][["names"]] <- tree_dendro$names
    cut_dendro <- cut_dendro(tree_dendro=tree_dendro, minClusterSize = minClusterSize, datan=datan, MEDissThres = MEDissThres,name_of_tree = paste0("Tree ",i,":"), plot = plot)
    res[[i]][["colors"]] <- cut_dendro$colors
    res[[i]][["MEs"]] <- cut_dendro$MEs
    i <- i+1
  }
  return(res)
}


#' Summarize results from a forest. Plot trees together.
#'
#' @param clust_res Clustering results from cut_trees call.
#' @return List
#' @export
nb_summary <- function(clust_res = NULL, plot = TRUE) {
  res <- list(
    dendros <- list(),
    names <- c(),
    colors <- c(),
    MEs <- data.frame(row.names = rownames(clust_res[[1]]$data))
  )
  
  n_MEs <- 0
  n_MEs_background <- 0
  for (tree in 1:length(clust_res)) {
    res$dendros[[tree]] <- clust_res[[tree]]$dendro
    res$names <- c(res$names, clust_res[[tree]]$names)
    tmp.col <- clust_res[[tree]]$colors
    tmp.col.new <- tmp.col
    j <- 1
    k <- -1
    for (col in unique(tmp.col)) {
      if (col != 0 | length(unique(tmp.col)) == 1) {
        tmp.col.new[tmp.col == col] <- n_MEs + j
        j <- j + 1
      }
      if (col == 0 & length(unique(tmp.col)) != 1) {
        tmp.col.new[tmp.col == col] <- -n_MEs_background + k
        k <- k - 1
      }
    }
    
    res$colors <- c(res$colors, tmp.col.new)
    if (ncol(clust_res[[tree]]$MEs) > 1) {
      tmp <- clust_res[[tree]]$MEs
      for (j in 1:ncol(tmp)) {
        if (colnames(tmp)[j] == "ME0") {
          n_MEs_background <- n_MEs_background + 1
          colnames(tmp)[j] <- paste0("ME0_", n_MEs_background)
        } else{
          n_MEs <- n_MEs +1
          colnames(tmp)[j] <- paste0("ME", (n_MEs))
        }
      }
      (dim(res$MEs))
      (dim(tmp))
      res$MEs <- cbind(res$MEs, tmp)
    } else{
      tmp <- clust_res[[tree]]$MEs
      colnames(tmp)[1] <- paste0("ME", (n_MEs + 1))
      res$MEs <- cbind(res$MEs, tmp)
      n_MEs <- n_MEs + 1
    }
  }
  
  cat("Netboost detected ",
      n_MEs,
      " modules and ",
      n_MEs_background,
      " background modules in ",
      length(clust_res),
      " trees.\n")
  cat("Average size of the modules was ", mean(table(res$colors[!(res$colors <=
                                                                    0)])), ".\n")
  cat(
    sum(res$colors <= 0),
    " of ",
    length(res$colors),
    " features (",
    (sum(res$colors <= 0) * 100 / length(res$colors)),
    "%) were not assigned to modules.\n"
  )
  
  if(plot == TRUE){
    colorHeight = 0.2
    layout(matrix(c(1:(
      2 * length(clust_res)
    )), nrow = 2), heights = c(1 - colorHeight, colorHeight))
  
    last_col <- 0
    plot_colors <- res$colors
    plot_colors[plots_colors <= 0] <- 0
    for (tree in 1:length(res$dendros)) {
      par(mar = c(0, 4, 8, 4))
      plot(res$dendro[[tree]],labels=FALSE)
      par(mar = c(4, 4, 0, 4))
      first_col <- last_col + 1
      last_col <- last_col + length(res$dendro[[tree]]$labels)
      plotColorUnderTree(res$dendro[[tree]], c(gray(level=0.7),rainbow(n = (length(unique(plot_colors))-1)))[plot_colors[first_col:last_col]+1])
    }
  }
  
  return(res)
}


#' Transfer of Netboost results to new data.
#' 
#' @param nb_summary Netboost results from nb_summary call.
#' @param new_data Data frame were rows correspond to samples and columns to features.
#' @return List
#' @export
nb_transfer <- function(nb_summary = NULL, new_data = NULL){
  if (!exists("new_data"))
    stop("datan must be provided")
  
  if (!(is.data.frame(new_data) && (nrow(new_data) > 0) && (ncol(new_data) > 0)))
    stop("new_data must be a data frame with dim() > (0,0).")
  
  if(length(nb_summary$colors) != ncol(new_data)){
    stop("The number of features in new_data must correspond to the number in nb_summary.")
  }
  
  if(sort(nb_summary$names) != sort(colnames(new_data))){
    stop("The features in new_data (colnames) must correspond to the features in nb_summary (nb_summary$names).")
  }
  
  new_data <- new_data[,nb_summary$names]
  # MEs <- data.frame(row.names = rownames(new_data), col.names = colnames(nb_summary$MEs))
  # # for(ME in colnames(nb_summary$MEs)){
  #   col <- if(strsplit(x = ME,split = "_")[1] == "ME0"){-strsplit(x = ME,split = "_")[2]}else{substring(x=ME, first = 3)}
  #   # MEs[,ME] <- new_data[,nb_summary$colors == col]
  #   
  # }
  background_MEs <- sum(unique(nb_summary$colors) <= 0)
  MEs <- moduleEigengenes(new_data, colors = nb_summary$colors+background_MEs)
  
  return(MEs)
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