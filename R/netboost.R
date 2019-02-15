# Load WGCNA, try to hide the welcome-message (only a try as it is
# printed...)
# Workaround using environment to force WGCNA skipping it's welcome-message.
Sys.setenv(ALLOW_WGCNA_THREADS=1)
suppressPackageStartupMessages(library(WGCNA))
Sys.unsetenv("ALLOW_WGCNA_THREADS")

library(colorspace)
library(parallel)

#' Netboost clustering.
#' 
#' The Netboost clustering is performed in three subsequent steps.
#' First, a filter of important edges in the network is calculated.
#' Next, pairwise distances are calculated.
#' Last, clustering is performed.
#' For details see Schlosser et al. doi...
#'
#' @param datan     Data frame were rows correspond to samples and columns to features.
#' @param stepno    Integer amount of boosting steps applied in the filtering step
#' @param until     Stop at index/column (if 0: iterate through all columns). For testing purposes in large datasets.
#' @param progress  Integer. If > 0, print progress after every X steps (mind: parallel!)
#' @param mode      Integer. Mode (0: x86, 1: FMA, 2: AVX). Features are only available
#'                  if compiled accordingly and available on the hardware.
#' @param softPower Integer. Exponent of the transformation. Set automatically based on the scale free topology criterion if unspecified.
#' @param max_singleton   Integer. The maximal singleton in the clustering. Usually equals the number of features.
#' @param plot      Logical. Should plots be created?
#' @param minClusterSize  Integer. The minimum number of features in one module.
#' @param MEDissThres Numeric. Module Eigengene Dissimilarity Threshold for merging close modules.
#' @param cores     Integer. Amount of CPU cores used (<=1 : sequential)
#' @param scale     Logical. Should data be scaled and centered?
#' @param verbose   Additional diagnostic messages.
#' @param nPC        Number of principal components and variance explained entries to be calculated. The number of returned variance explained entries is currently ‘min(nPC,10)’. If given ‘nPC’ is greater than 10, a warning is issued.
#' @param nb_min_varExpl        Minimum proportion of variance explained for returned module eigengenes. The number of PCs is capped at nPC.
#' @return dendros  A list of dendrograms. For each fully separate part of the network an individual dendrogram.
#' @return names    A vector of feature names.
#' @return colors   A vector of numeric color coding in matching order of names and module eigengene names (color = 3 -> variable in ME3).
#' @return MEs      Aggregated module measures (Module eigengenes).
#' @return varExplained  Proportion of variance explained per module eigengene.
#' @return dendros  A list of dendrograms. For each fully separate part of the network an individual dendrogram.
#' @export
netboost <- function(datan = NULL,
                     stepno = 20L, until = 0L,
                     progress = 1000L,
                     mode = 2L,
                     softPower = NULL,
                     max_singleton = dim(datan)[2],
                     plot = TRUE,
                     minClusterSize = 2L,
                     MEDissThres = 0.25,
                     nPC = 1,
                     nb_min_varExpl = 0.5,
                     cores = as.integer(getOption("mc.cores", 2)),
                     scale = TRUE,
                     verbose = getOption("verbose")) {
  # Initialize parallelization of WGCNA package.
  if (cores > 1) WGCNA::allowWGCNAThreads(nThreads = as.numeric(cores))

	if (ncol(datan)>5000000){
	   stop("A bug in sparse UPGMA currently prevents analyses with more than 5 million features.")
	}

	if(scale){
  message("Netboost: Scaling and centering data.")
  datan <- as.data.frame(scale(datan,center=TRUE,scale=TRUE))	
	}
  
  message("Netboost: Initialising filter step.")
  filter <- nb_filter(datan=datan, stepno=stepno, until=until, progress=progress, cores=cores,mode=mode)
  
  message("Netboost: Finished filter step.")
  
  if(is.null(softPower)){
    # Random subset out of allocation
    random_features <- sample(ncol(datan), min(c(10000,ncol(datan))))
    # Call the network topology analysis function
    sft <- pickSoftThreshold(datan[,random_features])
    softPower <- sft$powerEstimate
    message(paste0("Netboost: softPower was set to ", softPower, 
                   " based on the scale free topology criterion."))
  }
  
  message("Netboost: Initialising distance calculation.")
  dist <- nb_dist(datan=datan, filter=filter, softPower=softPower, cores=cores)
  message("Netboost: Finished distance calculation.")
  
  message("Netboost: Initialising clustering step.")
  results <- nb_clust(datan=datan, filter=filter, dist=dist, minClusterSize = minClusterSize, MEDissThres = MEDissThres, nPC = nPC, nb_min_varExpl = nb_min_varExpl,
                      max_singleton=max_singleton, cores=cores, plot = plot)
  message("Netboost: Finished clustering step.")
  
  message("Netboost: Finished Netboost.")
  return(results)
}

#' Calculate network adjacencies for filter
#' 
#' @param datan     Data frame were rows correspond to samples and columns to features.
#' @param filter    Filter-Matrix as generated by the nb_filter function.
#' @param softPower Integer. Exponent of the transformation. Set automatically based on the scale free topology criterion if unspecified.
#' @return Vector with adjacencies for the filter
calculate_adjacency <- function(datan=NULL,filter=NULL,softPower=2) {
  return(sapply(1:nrow(filter),
                function(i) {
                  abs(WGCNA::cor(datan[,filter[i,1]],datan[,filter[i,2]]))^softPower
                }))
}

#' Calculate distance
#' (external wrapper for internal C++ function)
#' Parallelisation inside C++ program with RcppParallel.
#' 
#' @param filter    Filter-Matrix as generated by the nb_filter function.
#' @param datan     Data frame were rows correspond to samples and columns to features.
#' @param softPower Integer. Exponent of the transformation. Set automatically based on the scale free topology criterion if unspecified.
#' @param cores     Integer. Amount of CPU cores used (<=1 : sequential).
#' @param verbose   Additional diagnostic messages.
#' @return Vector with distances.
#' @export
nb_dist <- function(filter=NULL,
                    datan=NULL,
                    softPower=2,
                    cores = getOption("mc.cores", 2L),
                    verbose = getOption("verbose")) {
  # if (is.null(filter) || is.null(adjacency))
  #   stop("Both filter and adjacency must be provided")

  if (!(is.matrix(filter) && (nrow(filter) > 0) && (ncol(filter) > 0)))
    stop("filter must be matrix with dim() > (0,0)")

  # if (!(is.vector(adjacency) && (length(adjacency) > 0)))
  #   stop("adjacency is required a vector with length > 0")
  cores <- max(cores, 1)

  ## RcppParallel amount of threads to be started
  setThreadOptions(numThreads = cores)

  return(cpp_dist_tom(filter,
                      calculate_adjacency(datan=datan,
                                          filter=filter,
                                          softPower=softPower)))
}

#' Calculate dendrogram for a sparse distance matrix
#' (external wrapper MC-UPGMA clustering package Loewenstein et al.
#' 
#' @param filter    Filter-Matrix as generated by the nb_filter function.
#' @param dist      Distance-Matrix as generated by the nb_dist function.
#' @param max_singleton   Integer. The maximal singleton in the clustering. Usually equals the number of features.
#' @param cores     Integer. Amount of CPU cores used (<=1 : sequential)
#' @param verbose   Logical. Additional diagnostic messages.
#' @return Raw dendrogram to be processed by tree_search and tree_dendro.
nb_mcupgma <- function(filter = NULL,
                       dist = NULL,
                       max_singleton = NULL,
                       cores = getOption("mc.cores", 2L),
                       verbose = getOption("verbose")) {
  # Deletes all files under netboostTmpPath(), esp. clustering/iteration_
  netboostTmpCleanup()
  
  	if (max_singleton>5000000){
	   stop("A bug in sparse UPGMA currently prevents analyses with more than 5 million features.")
	}

  if (!dir.create(file.path(netboostTmpPath(), "clustering")))
    stop(paste("Unable to create:", file.path(netboostTmpPath(), "clustering")))
  
  file_dist_edges <- file.path(netboostTmpPath(), "clustering", "dist.edges")
  file_dist_tree <- file.path(netboostTmpPath(), "clustering", "dist.mcupgma_tree")

#  write.table(file="clustering/dist.edges",
  write.table(file=file_dist_edges,
              cbind(format(filter, scientific=FALSE, trim=TRUE),
                    format(dist, scientific=FALSE, trim=TRUE)),
              row.names=FALSE,
              col.names=FALSE,
              sep="\t",
              quote=FALSE)
  
#  system("gzip -f clustering/dist.edges")
  ret <- system2("gzip",
                 args = c("-f", file_dist_edges,
                          ifelse(verbose, "--verbose", "")))

  ## If gzip compressed file, the original file (and variable) is replaced
  if (ret == 0 && file.exists(paste0(file_dist_edges, ".gz")))
    file_dist_edges <- paste0(file_dist_edges, ".gz")
  else
    warning(paste("Gzip maybe failed on:", file_dist_edges, "Return:", ret))
  
  ret <- mcupgma_exec(exec="cluster.pl",
                      "-max_distance", 1,
                      "-max_singleton", max_singleton,
                      "-iterations 1000 -heap_size 10000000 -num_hash_buckets 40",
                      "-jobs", cores,
                      "-retries 1",
                      "-output_tree_file", file_dist_tree,
                      "-split_unmodified_edges",
                      cores,     # Ist das hier richtig?
                      file_dist_edges,
                      console = FALSE)
  
  if (verbose)
    print(ret)
  
  if (!file.exists(file_dist_tree) || file.info(file_dist_tree)$size == 0)
    stop("No output file created. mcupgma error :(")

  return(as.matrix(read.table(file=file_dist_tree,
                              row.names=NULL,
                              col.names=c("cluster_id1","cluster_id2","distance","cluster_id3"))))
}

#' Extracts independent trees from nb_mcupgma results
#' (external wrapper for internal C++ function)
#'
#' @param forest Raw dendrogram-matrix as generated by the nb_mcupgma function.
#' @return List
tree_search <- function(forest=NULL) {
  # Check for integer values cannot be done here, as either the user must
  # have set up all values with as.integer() or R delivers default numeric
  # (double). In that case, Rcpp converts the matrix.
  # (Alternative: convert manually "matrix(as.integer(forest), nrow=nrow(forest))")
  forest <- as.matrix(forest)
  
  if (is.null(forest) || !is.matrix(forest))
    stop("forest must be provided (as integer matrix)")
  
  return(cpp_tree_search(forest))
}


#' Calculate the dendrogram for an individual tree
#'
#' @param tree A list with two elements. ids, which is an integer vector of feature identifiers and rows, which is an integer vector of selected rows in the corresponding forest
#' @param datan     Data frame were rows correspond to samples and columns to features.
#' @param forest Raw dendrogram-matrix as generated by the nb_mcupgma function.
#' @return List of tree specific objects including dendrogram, tree data and features.
tree_dendro <- function(tree=NULL, datan=NULL, forest=NULL) {
  index.features <- tree$ids[tree$ids <= dim(datan)[2]]
  data_tree <- datan[,index.features]
  
  colnames_tree <- colnames(datan)[index.features]
  tree_cluster <- forest[tree$rows,,drop=FALSE]

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
  dendro$merge <- tree_cluster[,c(1,2),drop=FALSE]
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
  
#' Module detection for an individual tree
#'
#' @param tree_dendro List of tree specific objects including dendrogram, tree data and features originating from the tree_dendro function.
#' @param minClusterSize  Integer. The minimum number of features in one module.
#' @param datan     Data frame were rows correspond to samples and columns to features.
#' @param MEDissThres Numeric. Module Eigengene Dissimilarity Threshold for merging close modules.
#' @param name_of_tree String. Annotating plots and messages.
#' @param plot      Logical. Should plots be created?
#' @param nPC        Number of principal components and variance explained entries to be calculated. The number of returned variance explained entries is currently ‘min(nPC,10)’. If given ‘nPC’ is greater than 10, a warning is issued.
#' @param nb_min_varExpl        Minimum proportion of variance explained for returned module eigengenes. The number of PCs is capped at nPC.
#' @return Object of class hclust
cut_dendro <- function(tree_dendro=NULL, minClusterSize= 2L, 
                       datan=NULL, MEDissThres = NULL,
                       name_of_tree="", plot = TRUE, nPC = 1, nb_min_varExpl = 0.5) {
  dynamicMods <- cutreeDynamic(dendro = tree_dendro$dendro, method="tree", deepSplit = TRUE, minClusterSize = minClusterSize)
  ### Merging of Dynamic Modules ###
  # Calculate eigengenes
  MEList <- nb_moduleEigengenes(expr=tree_dendro$data, colors = dynamicMods, nPC = nPC, nb_min_varExpl = nb_min_varExpl)
  MEs <- MEList$nb_eigengenes
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
    MEList <- nb_moduleEigengenes(expr=tree_dendro$data, colors = merged$colors, nPC = nPC, nb_min_varExpl = nb_min_varExpl)
    MEs <- MEList$nb_eigengenes
    MEDiss <- 1-cor(MEs);
    if(length(MEDiss) > 1){
      METree <- hclust(as.dist(MEDiss), method = "average");
      if(plot == TRUE & length(tree_dendro$dendro$labels)>2){
        plot(METree, main =paste0(name_of_tree,"Clustering of merged module eigengenes"),
           xlab = "", sub = "")
        plotDendroAndColors(dendro=tree_dendro$dendro, colors=mergedColors,"Merged Dynamic", dendroLabels = FALSE, hang = 0.01, addGuide = TRUE, guideHang = 0.05, main=paste0(name_of_tree,"Cluster Dendrogram"))
      }
    }
    }else{
    cat("\nOnly one module in ",name_of_tree,".\n")
    mergedColors <- dynamicMods
    if(length(tree_dendro$dendro$labels)>2){
      if(plot == TRUE){
        plot(tree_dendro$dendro, main=paste0(name_of_tree,"Cluster Dendrogram (Tree maximaly consists out of one module.)"))
      }
    }else{cat("\nOnly two elements in the one module in ",name_of_tree," (no plot generated).\n")}
  }
  cat("\nNetboost extracted",length(table(mergedColors)),"modules (including background) with an average size of",mean(table(mergedColors)[-1])," (excluding background) from ",substr(name_of_tree,start=1,stop=(nchar(name_of_tree)-1)),".\n")
  return(list(colors=mergedColors,MEs=MEs,varExplained=MEList$varExplained,nb_rotations=MEList$nb_rotations))
}

#' Module detection for the results from a nb_mcupgma call
#'
#' @param trees List of trees, where one tree is a list of ids and rows
#' @param datan     Data frame were rows correspond to samples and columns to features.
#' @param forest Raw dendrogram-matrix as generated by the nb_mcupgma function.
#' @param minClusterSize  Integer. The minimum number of features in one module.
#' @param MEDissThres Numeric. Module Eigengene Dissimilarity Threshold for merging close modules.
#' @param plot      Logical. Should plots be created?
#' @param nPC        Number of principal components and variance explained entries to be calculated. The number of returned variance explained entries is currently ‘min(nPC,10)’. If given ‘nPC’ is greater than 10, a warning is issued.
#' @param nb_min_varExpl        Minimum proportion of variance explained for returned module eigengenes. The number of PCs is capped at nPC.
#' @return List
cut_trees <- function(trees=NULL,
                      datan=NULL, 
                      forest=NULL,
                      minClusterSize= 2L,
                      MEDissThres = NULL,
                      plot = TRUE,
                      nPC = 1,
                      nb_min_varExpl = 0.5) {
  res <- list()
  i <- 1L
  for(tree in trees){
    tree_dendro <- netboost:::tree_dendro(tree=tree,datan=datan,forest=forest)
    res[[i]] <- list()
    res[[i]][["dendro"]] <- tree_dendro$dendro
    res[[i]][["data"]] <- tree_dendro$data
    res[[i]][["names"]] <- tree_dendro$names
    cut_dendro <- netboost:::cut_dendro(tree_dendro=tree_dendro, minClusterSize = minClusterSize, datan=datan, MEDissThres = MEDissThres,name_of_tree = paste0("Tree ",i,":"), plot = plot, nPC = nPC, nb_min_varExpl = nb_min_varExpl)
    res[[i]][["colors"]] <- cut_dendro$colors
    res[[i]][["MEs"]] <- cut_dendro$MEs
    res[[i]][["varExplained"]] <- cut_dendro$varExplained
    res[[i]][["nb_rotations"]] <- cut_dendro$nb_rotations
    i <- i+1
  }
  return(res)
}

#' Netboost clustering step
#'
#' @param filter    Filter-Matrix as generated by the nb_filter function.
#' @param dist      Distance-Matrix as generated by the nb_dist function.
#' @param datan     Data frame were rows correspond to samples and columns to features.
#' @param max_singleton   Integer. The maximal singleton in the clustering. Usually equals the number of features.
#' @param minClusterSize  Integer. The minimum number of features in one module.
#' @param MEDissThres Numeric. Module Eigengene Dissimilarity Threshold for merging close modules.
#' @param cores     Integer. Amount of CPU cores used (<=1 : sequential)
#' @param plot Logical. Create plot.
#' @param nPC        Number of principal components and variance explained entries to be calculated. The number of returned variance explained entries is currently ‘min(nPC,10)’. If given ‘nPC’ is greater than 10, a warning is issued.
#' @param nb_min_varExpl        Minimum proportion of variance explained for returned module eigengenes. The number of PCs is capped at nPC.
#' @return List
#' @export
nb_clust <- function(filter = NULL,
                     dist = NULL,
                     datan = NULL,
                     max_singleton = dim(datan)[2],
                     minClusterSize = 2L,
                     MEDissThres = 0.25,
                     cores = getOption("mc.cores", 2L),
                     plot = TRUE,
                     nPC = 1,
                     nb_min_varExpl = 0.5) {
  forest <- netboost:::nb_mcupgma(filter=filter,dist=dist,max_singleton=max_singleton,cores=cores)
  trees <- netboost:::tree_search(forest)
  results <- netboost:::cut_trees(trees=trees,datan=datan, forest=forest, minClusterSize = minClusterSize, MEDissThres = MEDissThres, plot = plot, nPC = nPC, nb_min_varExpl = nb_min_varExpl)
  sum_res <- netboost:::nb_summary(clust_res = results, plot = plot)
  return(sum_res)
}

#' Summarize results from a forest. Plot trees together.
#'
#' @param clust_res Clustering results from cut_trees call.
#' @param plot Logical. Create plot.
#' @return List
nb_summary <- function(clust_res = NULL, plot = TRUE) {
  res <- vector("list")
  n_MEs <- 0
  n_MEs_background <- 0
  for (tree in 1:length(clust_res)) {
    res$dendros[[tree]] <- clust_res[[tree]]$dendro
    res$names <- c(res$names, clust_res[[tree]]$names)
    tmp.col <- clust_res[[tree]]$colors
    tmp.col.new <- tmp.col
    tmp_MEs <- clust_res[[tree]]$MEs
    tmp_MEs_new <- tmp_MEs
    for (col in unique(tmp.col)) {
      if (col != 0 | length(unique(tmp.col)) == 1) {
        n_MEs <- n_MEs +1
        tmp.col.new[tmp.col == col] <- n_MEs
        colnames(tmp_MEs_new)[grepl(pattern=paste0("ME",col,"_"),colnames(tmp_MEs))] <- gsub(pattern=paste0("ME",col,"_"),replacement=paste0("ME",(n_MEs),"_"),colnames(tmp_MEs_new)[grepl(pattern=paste0("ME",col,"_"),colnames(tmp_MEs))])
      }
      if (col == 0 & length(unique(tmp.col)) != 1) {
        n_MEs_background <- n_MEs_background + 1
        tmp.col.new[tmp.col == col] <- -n_MEs_background
        colnames(tmp_MEs_new)[grepl(pattern="ME0_",colnames(tmp_MEs))] <- gsub(pattern="ME0_",replacement=paste0("ME0_", n_MEs_background,"_"),colnames(tmp_MEs_new)[grepl(pattern="ME0_",colnames(tmp_MEs))])
      }
    }
    res$colors <- c(res$colors, tmp.col.new)
    if("MEs" %in% names(res)){res$MEs <- cbind(res$MEs, tmp_MEs_new)}else{res$MEs <- tmp_MEs_new}
    if("varExplained" %in% names(res)){
      res$varExplained <- cbind(res$varExplained , clust_res[[tree]]$varExplained)
    }else{
        res$varExplained <-  clust_res[[tree]]$varExplained
    }
  }
  rownames(res$varExplained)<- paste0("PC",1:nrow(res$varExplained))
  colnames(res$varExplained)<- unique(unlist(lapply(strsplit(split="_pc",colnames(res$MEs)),FUN=function(x){x[1]})))

  cat("Netboost detected ",
      n_MEs,
      " modules and ",
      n_MEs_background,
      " background modules in ",
      length(clust_res),
      " trees.\n")
  cat("Average size of the modules was ", mean(table(res$colors[!(res$colors <= 0)])), ".\n")
  cat(
    sum(res$colors <= 0),
    " of ",
    length(res$colors),
    " features (",
    (sum(res$colors <= 0) * 100 / length(res$colors)),
    "%) were not assigned to modules.\n"
  )
  
  if(plot == TRUE){
    nb_plot_dendro(nb_summary = res,labels=FALSE)
  }
  return(res)
}


#' Transfer of Netboost clustering to new data.
#' 
#' @param nb_summary Netboost results as generated by the nb_summary function.
#' @param new_data Data frame were rows correspond to samples and columns to features.
#' @param scale     Logical. Should data be scaled and centered?
#' @return List
#' @export
nb_transfer <- function(nb_summary = NULL, new_data = NULL, scale = FALSE){
  if (!exists("new_data"))
    stop("datan must be provided")
  
  if (!(is.data.frame(new_data) && (nrow(new_data) > 0) && (ncol(new_data) > 0)))
    stop("new_data must be a data frame with dim() > (0,0).")
  
  if(length(nb_summary$colors) != ncol(new_data)){
    stop("The number of features in new_data must correspond to the number in nb_summary.")
  }
  
  if(!identical(sort(nb_summary$names), sort(colnames(new_data)))){
    stop("The features in new_data (colnames) must correspond to the features in nb_summary (nb_summary$names).")
  }
  
  new_data <- new_data[,nb_summary$names]
  
  if(scale){
	new_data <- as.data.frame(scale(new_data,center=TRUE,scale=TRUE))
  }
  
  #will replace with rotation matrix transfer
  MEs <- nb_moduleEigengenes(expr=new_data, colors = nb_summary$colors)$nb_eigengenes

  colnames(MEs)[lapply(strsplit(x = colnames(MEs), split = "-"),FUN = length) > 1] <- paste0("ME0_",substring(text = colnames(MEs)[lapply(strsplit(x = colnames(MEs), split = "-"),FUN = length) > 1], first = 4))
  MEs <- MEs[,colnames(nb_summary$MEs)]
  rownames(MEs) <- rownames(new_data)
  return(MEs)
}


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
  
  if (ncol(datan)>5000000){
    stop("A bug in sparse UPGMA currently prevents analyses with more than 5 million features.")
  }
  
  message(paste("Netboost: Filtering"))
  
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

#' TCGA RNA and methylation measurement on chromosome 18 for 180 AML patients.
#'
#' @format A data frame with 180 rows and 5283 variables:
#' @source \url{http://www.tcga.com/}
"tcga_aml_meth_rna_chr18"


#' Plot dendrogram from Netboost output.
#' 
#' @param nb_summary Netboost results as generated by the nb_summary function.
#' @param labels TRUE/FALSE indicator of whether labels should be attached to the leafs.
#' @param main Plot title.
#' @param colorsrandom TRUE/FALSE indicator of whether module colors should be shuffeled.
#' @param seed Seed for shuffeling of module colors.
#' @export
nb_plot_dendro <- function(nb_summary = NULL,labels=FALSE,main="",colorsrandom=FALSE,seed=NULL){
  if (!exists("nb_summary"))
    stop("Netboost output (nb_summary) must be provided.")
  
  colorHeight = 0.2
  layout(matrix(c(1:(
    2 * length(nb_summary$dendros)
  )), nrow = 2), heights = c(1 - colorHeight, colorHeight))
  
  last_col <- 0
  n_colors <- length(unique(nb_summary$colors))
  middle <- floor(n_colors/2)
  if(colorsrandom){
    if(!is.null(seed)){
      set.seed(seed)
    }
    shuffel_index <- sample(x = n_colors,size = n_colors)
  }else{
    shuffel_index <- 1:n_colors
  }
  shuffel_index <- c(shuffel_index[(2*middle+1)][(2*middle+1)==n_colors],rbind(shuffel_index[1:middle],shuffel_index[(middle+1):(2*middle)]))
  plot_colors <- colorspace::rainbow_hcl(n = (length(unique(nb_summary$colors))))[shuffel_index][as.factor(nb_summary$colors)]
  plot_colors[nb_summary$colors <= 0] <- gray(level = 0.7)
  for (tree in 1:length(nb_summary$dendros)) {
    par(mar = c(0, 4, 8, 4))
    first_col <- last_col + 1
    last_col <- last_col + length(nb_summary$dendros[[tree]]$labels)
    if(labels){
      plot(nb_summary$dendros[[tree]],labels=nb_summary$names[first_col:last_col],main=main)
    }else{
      plot(nb_summary$dendros[[tree]],labels=FALSE,main=main)
    }
    par(mar = c(4, 4, 0, 4))
    WGCNA::plotColorUnderTree(nb_summary$dendros[[tree]], colors=plot_colors[first_col:last_col],rowLabels="")
  }
}


#' Netboost module aggregate extraction.
#' 
#' This is a modification of WGCNA:::moduleEigengenes() (version WGCNA_1.66) to include more than the first principal component.
#' For details see WGCNA:::moduleEigengenes().
#'
#' @param expr     Expression data for a single set in the form of a data frame where rows are samples and columns are genes (probes).
#' @param colors    A vector of the same length as the number of probes in ‘expr’, giving module color for all probes (genes). Color ‘"grey"’ is reserved for unassigned genes.     Expression
#' @param impute    If ‘TRUE’, expression data will be checked for the presence of ‘NA’ entries and if the latter are present, numerical data will be imputed, using function ‘impute.knn’ and probes from the same module as the missing datum. The function ‘impute.knn’ uses a fixed random seed giving repeatable results.
#' @param nPC       Number of principal components and variance explained entries to be calculated. The number of returned variance explained entries is currently ‘min(nPC,10)’. If given ‘nPC’ is greater than 10, a warning is issued.
#' @param align     Controls whether eigengenes, whose orientation is undetermined, should be aligned with average expression (‘align = "along average"’, the default) or left as they are (‘align = ""’). Any other value will trigger an error.
#' @param excludeGrey   Should the improper module consisting of 'grey' genes be excluded from the eigengenes?
#' @param grey          Value of ‘colors’ designating the improper module. Note that if ‘colors’ is a factor of numbers, the default value will be incorrect.
#' @param subHubs       Controls whether hub genes should be substituted for missing eigengenes. If ‘TRUE’, each missing eigengene (i.e., eigengene whose calculation failed and the error was trapped) will be replaced by a weighted average of the most connected hub genes in the corresponding module. If this calculation fails, or if ‘subHubs==FALSE’, the value of ‘trapErrors’ will determine whether the offending module will be removed or whether the function will issue an error and stop.
#' @param trapErrors    Controls handling of errors from that may arise when there are too many ‘NA’ entries in expression data. If ‘TRUE’, errors from calling these functions will be trapped without abnormal exit.  If ‘FALSE’, errors will cause the function to stop. Note, however, that ‘subHubs’ takes precedence in the sense that if ‘subHubs==TRUE’ and ‘trapErrors==FALSE’, an error will be issued only if both the principal component and the hubgene calculations have failed.
#' @param returnValidOnly   logical; controls whether the returned data frame of module eigengenes contains columns corresponding only to modules whose eigengenes or hub genes could be calculated correctly (‘TRUE’), or whether the data frame should have columns for each of the input color labels (‘FALSE’).
#' @param softPower     The power used in soft-thresholding the adjacency matrix. Only used when the hubgene approximation is necessary because the principal component calculation failed. It must be non-negative. The default value should only be changed if there is a clear indication that it leads to incorrect results.
#' @param scale         logical; can be used to turn off scaling of the expression data before calculating the singular value decomposition. The scaling should only be turned off if the data has been scaled previously, in which case the function can run a bit faster. Note however that the function first imputes, then scales the expression data in each module. If the expression contain missing data, scaling outside of the function and letting the function impute missing data may lead to slightly different results than if the data is scaled within the function.
#' @param verbose       Controls verbosity of printed progress messages. 0 means silent, up to (about) 5 the verbosity gradually increases.
#' @param indent        A single non-negative integer controlling indentation of printed messages. 0 means no indentation, each unit above that adds two spaces.
#' @param nb_min_varExpl        Minimum proportion of variance explained for returned module eigengenes. Is capped at nPC.
#' 
#' @return eigengenes   Module eigengenes in a dataframe, with each column corresponding to one eigengene. The columns are named by the corresponding color with an ‘"ME"’ prepended, e.g., ‘MEturquoise’ etc. If ‘returnValidOnly==FALSE’, module eigengenes whose calculation failed have all components set to ‘NA’.
#' @return averageExpr  If ‘align == "along average"’, a dataframe containing average normalized expression in each module. The columns are named by the corresponding color with an ‘"AE"’ prepended, e.g., ‘AEturquoise’ etc.
#' @return varExplained A dataframe in which each column corresponds to a module, with the component ‘varExplained[PC, module]’ giving the variance of module ‘module’ explained by the principal component no. ‘PC’. The calculation is exact irrespective of the number of computed principal components. At most 10 variance explained values are recorded in this dataframe.
#' @return nPC          A copy of the input ‘nPC’.
#' @return validMEs     A boolean vector. Each component (corresponding to the columns in ‘data’) is ‘TRUE’ if the corresponding eigengene is valid, and ‘FALSE’ if it is invalid. Valid eigengenes include both principal components and their hubgene approximations. When ‘returnValidOnly==FALSE’, by definition all returned eigengenes are valid and the entries of ‘validMEs’ are all ‘TRUE’.
#' @return validColors  A copy of the input colors with entries corresponding to invalid modules set to ‘grey’ if given, otherwise 0 if ‘colors’ is numeric and "grey" otherwise.
#' @return allOK        Boolean flag signalling whether all eigengenes have been calculated correctly, either as principal components or as the hubgene average approximation.
#' @return allPC        Boolean flag signalling whether all returned eigengenes are principal components.
#' @return isPC         Boolean vector. Each component (corresponding to the columns in ‘eigengenes’) is ‘TRUE’ if the corresponding eigengene is the first principal component and ‘FALSE’ if it is the hubgene approximation or is invalid.
#' @return isHub        Boolean vector. Each component (corresponding to the columns in ‘eigengenes’) is ‘TRUE’ if the corresponding eigengene is the hubgene approximation and ‘FALSE’ if it is the first principal component or is invalid.
#' @return validAEs     Boolean vector. Each component (corresponding to the columns in ‘eigengenes’) is ‘TRUE’ if the corresponding module average expression is valid.
#' @return allAEOK      Boolean flag signalling whether all returned module average expressions contain valid data. Note that ‘returnValidOnly==TRUE’ does not imply ‘allAEOK==TRUE’: some invalid average expressions may be returned if their corresponding eigengenes have been calculated correctly.
#' @export
nb_moduleEigengenes <- function (expr, colors, impute = TRUE, nPC = 1, align = "along average", 
    excludeGrey = FALSE, grey = if (is.numeric(colors)) 0 else "grey", 
    subHubs = TRUE, trapErrors = FALSE, returnValidOnly = trapErrors, 
    softPower = 6, scale = TRUE, verbose = 0, indent = 0,nb_min_varExpl=0.5) 
{
    spaces = indentSpaces(indent)
    if (verbose == 1) 
        printFlush(paste(spaces, "moduleEigengenes: Calculating", 
            nlevels(as.factor(colors)), "module eigengenes in given set."))
    if (is.null(expr)) {
        stop("moduleEigengenes: Error: expr is NULL. ")
    }
    if (is.null(colors)) {
        stop("moduleEigengenes: Error: colors is NULL. ")
    }
    if (is.null(dim(expr)) || length(dim(expr)) != 2) 
        stop("moduleEigengenes: Error: expr must be two-dimensional.")
    if (dim(expr)[2] != length(colors)) 
        stop("moduleEigengenes: Error: ncol(expr) and length(colors) must be equal (one color per gene).")
    if (is.factor(colors)) {
        nl = nlevels(colors)
        nlDrop = nlevels(colors[, drop = TRUE])
        if (nl > nlDrop) 
            stop(paste("Argument 'colors' contains unused levels (empty modules). ", 
                "Use colors[, drop=TRUE] to get rid of them."))
    }
    if (softPower < 0) 
        stop("softPower must be non-negative")
    alignRecognizedValues = c("", "along average")
    if (!is.element(align, alignRecognizedValues)) {
        printFlush(paste("ModulePrincipalComponents: Error:", 
            "parameter align has an unrecognised value:", align, 
            "; Recognized values are ", alignRecognizedValues))
        stop()
    }
    maxVarExplained = 10 
    if (nPC > maxVarExplained) 
        warning(paste("Given nPC is too large. Will use value", 
            maxVarExplained))
    nVarExplained = min(nPC, maxVarExplained)
    modlevels = levels(factor(colors))
    if (excludeGrey) 
        if (sum(as.character(modlevels) != as.character(grey)) > 
            0) {
            modlevels = modlevels[as.character(modlevels) != 
                as.character(grey)]
        }
        else {
            stop(paste("Color levels are empty. Possible reason: the only color is grey", 
                "and grey module is excluded from the calculation."))
        }
    PrinComps = data.frame(matrix(NA, nrow = dim(expr)[[1]], 
        ncol = length(modlevels)))
    nb_PrinComps <- data.frame(matrix(NA, nrow = dim(expr)[[1]], 
        ncol = 0))
    nb_rotations <- vector("list", length(modlevels))
    averExpr = data.frame(matrix(NA, nrow = dim(expr)[[1]], ncol = length(modlevels)))
    varExpl = data.frame(matrix(NA, nrow = nVarExplained, ncol = length(modlevels)))
    validMEs = rep(TRUE, length(modlevels))
    validAEs = rep(FALSE, length(modlevels))
    isPC = rep(TRUE, length(modlevels))
    isHub = rep(FALSE, length(modlevels))
    validColors = colors
    names(PrinComps) = paste(moduleColor.getMEprefix(), modlevels, 
        sep = "")
    names(averExpr) = paste("AE", modlevels, sep = "")
    for (i in c(1:length(modlevels))) {
        if (verbose > 1) 
            printFlush(paste(spaces, "moduleEigengenes : Working on ME for module", 
                modlevels[i]))
        modulename = modlevels[i]
        restrict1 = as.character(colors) == as.character(modulename)
        if (verbose > 2) 
            printFlush(paste(spaces, " ...", sum(restrict1), 
                "genes"))
        datModule = as.matrix(t(expr[, restrict1]))
        n = dim(datModule)[1]
        p = dim(datModule)[2]
        pc = try({
            if (nrow(datModule) > 1 && impute) {
                seedSaved = FALSE
                if (exists(".Random.seed")) {
                  saved.seed = .Random.seed
                  seedSaved = TRUE
                }
                if (any(is.na(datModule))) {
                  if (verbose > 5) 
                    printFlush(paste(spaces, " ...imputing missing data"))
                  datModule = impute.knn(datModule, k = min(10, 
                    nrow(datModule) - 1))
                  try({
                    if (!is.null(datModule$data)) 
                      datModule = datModule$data
                  }, silent = TRUE)
                }
                if (seedSaved) 
                  .Random.seed <<- saved.seed
            }
            if (verbose > 5) 
                printFlush(paste(spaces, " ...scaling"))
            if (scale) 
                datModule = t(scale(t(datModule)))
            if (verbose > 5) 
                printFlush(paste(spaces, " ...calculating SVD"))
            svd1 = svd(datModule, nu = min(n, p, nPC), nv = min(n, 
                p, nPC))
            nb_PCA <- prcomp(x=t(datModule), retx = TRUE, center = FALSE, scale. = FALSE,tol = NULL, rank. = NULL)
#            print(dim(nb_PCA$x))
#            nb_PCs <- nb_PCA$x
#            nb_rotation <- nb_PCA$rotation
            if (verbose > 5) 
                printFlush(paste(spaces, " ...calculating PVE"))
            veMat = cor(svd1$v[, c(1:min(n, p, nVarExplained))], 
                t(datModule), use = "p")
            varExpl[c(1:min(n, p, nVarExplained)), i] = rowMeans(veMat^2, 
                na.rm = TRUE)
            svd1$v[, 1]
        }, silent = TRUE)
        if (class(pc) == "try-error") {
            if ((!subHubs) && (!trapErrors)) 
                stop(pc)
            if (subHubs) {
                if (verbose > 0) {
                  printFlush(paste(spaces, " ..principal component calculation for module", 
                    modulename, "failed with the following error:"))
                  printFlush(paste(spaces, "     ", pc, spaces, 
                    " ..hub genes will be used instead of principal components."))
                }
                isPC[i] = FALSE
                pc = try({
                  scaledExpr = scale(t(datModule))
                  covEx = cov(scaledExpr, use = "p")
                  covEx[!is.finite(covEx)] = 0
                  modAdj = abs(covEx)^softPower
                  kIM = (rowMeans(modAdj, na.rm = TRUE))^3
                  if (max(kIM, na.rm = TRUE) > 1) 
                    kIM = kIM - 1
                  kIM[is.na(kIM)] = 0
                  hub = which.max(kIM)
                  alignSign = sign(covEx[, hub])
                  alignSign[is.na(alignSign)] = 0
                  isHub[i] = TRUE
                  pcxMat = scaledExpr * matrix(kIM * alignSign, 
                    nrow = nrow(scaledExpr), ncol = ncol(scaledExpr), 
                    byrow = TRUE)/sum(kIM)
                  pcx = rowMeans(pcxMat, na.rm = TRUE)
                  varExpl[1, i] = mean(cor(pcx, t(datModule), 
                    use = "p")^2, na.rm = TRUE)
                  pcx
                }, silent = TRUE)
            }
        }
        if (class(pc) == "try-error") {
            if (!trapErrors) 
                stop(pc)
            if (verbose > 0) {
                printFlush(paste(spaces, " ..ME calculation of module", 
                  modulename, "failed with the following error:"))
                printFlush(paste(spaces, "     ", pc, spaces, 
                  " ..the offending module has been removed."))
            }
            warning(paste("Eigengene calculation of module", 
                modulename, "failed with the following error \n     ", 
                pc, "The offending module has been removed.\n"))
            validMEs[i] = FALSE
            isPC[i] = FALSE
            isHub[i] = FALSE
            validColors[restrict1] = grey
        }
        else {
            PrinComps[, i] = pc
            nb_nPCs <- min(c(which(cumsum(varExpl[c(1:min(n, p, nVarExplained)), i])>nb_min_varExpl),nVarExplained))
            nb_PrinComps <- cbind(nb_PrinComps,nb_PCA$x[,1:nb_nPCs])
            colnames(nb_PrinComps)[(ncol(nb_PrinComps)-nb_nPCs+1):ncol(nb_PrinComps)] <- paste0(moduleColor.getMEprefix(), modlevels[i],"_pc",1:nb_nPCs)
            nb_rotations[[i]] <- nb_PCA$rotations[,1:nb_nPCs]
            ae = try({
                if (isPC[i]) 
                  scaledExpr = scale(t(datModule))
                averExpr[, i] = rowMeans(scaledExpr, na.rm = TRUE)
                if (align == "along average") {
                  if (verbose > 4) 
                    printFlush(paste(spaces, " .. aligning module eigengene with average expression."))
                  corAve = cor(averExpr[, i], PrinComps[, i], 
                    use = "p")
                  if (!is.finite(corAve)) 
                    corAve = 0
                  if (corAve < 0) 
                    PrinComps[, i] = -PrinComps[, i]
                }
                0
            }, silent = TRUE)
            if (class(ae) == "try-error") {
                if (!trapErrors) 
                  stop(ae)
                if (verbose > 0) {
                  printFlush(paste(spaces, " ..Average expression calculation of module", 
                    modulename, "failed with the following error:"))
                  printFlush(paste(spaces, "     ", ae, spaces, 
                    " ..the returned average expression vector will be invalid."))
                }
                warning(paste("Average expression calculation of module", 
                  modulename, "failed with the following error \n     ", 
                  ae, "The returned average expression vector will be invalid.\n"))
            }
            validAEs[i] = !(class(ae) == "try-error")
        }
    }
    allOK = (sum(!validMEs) == 0)
    if (returnValidOnly && sum(!validMEs) > 0) {
        PrinComps = PrinComps[, validMEs]
        averExpr = averExpr[, validMEs]
        varExpl = varExpl[, validMEs]
        validMEs = rep(TRUE, times = ncol(PrinComps))
        isPC = isPC[validMEs]
        isHub = isHub[validMEs]
        validAEs = validAEs[validMEs]
    }
    allPC = (sum(!isPC) == 0)
    allAEOK = (sum(!validAEs) == 0)
    list(eigengenes = PrinComps, averageExpr = averExpr, varExplained = varExpl, 
        nPC = nPC, validMEs = validMEs, validColors = validColors, 
        allOK = allOK, allPC = allPC, isPC = isPC, isHub = isHub, 
        validAEs = validAEs, allAEOK = allAEOK,nb_eigengenes=nb_PrinComps,nb_rotations=nb_rotations)
}
