# Load WGCNA, try to hide the annoying mega-message (only a try as it is
# printed...)
# Workaround using environment to force WGCNA skipping it's mega-message.
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
#' @param softPower Integer. Exponent of the transformation
#' @param cores Integer. Amount of CPU cores used (<=1 : sequential)
#' @param datan     Data frame were rows correspond to samples and columns to features.
#' @param stepno    Integer amount of boosting steps applied in the filtering step
#' @param until     Stop at index/column (if 0: iterate through all columns)
#' @param progress  Integer. If > 0, print progress after every X steps (mind: parallel!)
#' @param mode      Integer. Mode (0: x86, 1: FMA, 2: AVX). Features are only available
#'                  if compiled accordingly and available on the hardware.
#' @param max_singleton   Integer. The maximal singleton in the clustering. Usually equals the number of features.
#' @param minClusterSize  Integer. The minimum number of features in one module.
#' @param MEDissThres Numeric. Module Eigengene Dissimilarity Threshold for merging close modules.
#' @param plot Logical. Create plot.
#' @param verbose   Additional diagnostic messages.
#' @return List
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
                     cores = as.integer(getOption("mc.cores", 2)),
                     verbose = getOption("verbose")) {
  # Initialize parallelization of WGCNA package.
  if (cores > 1) WGCNA::allowWGCNAThreads(nThreads = as.numeric(cores))

	if (ncol(datan)>5000000){
	   stop("A bug in sparse UPGMA currently prevents analyses with more than 5 million features.")
	}

  message("Netboost: Scaling and centering data.")
  datan <- as.data.frame(scale(datan,center=TRUE,scale=TRUE))
  
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
  results <- nb_clust(datan=datan, filter=filter, dist=dist, minClusterSize = minClusterSize, MEDissThres = MEDissThres,
                      max_singleton=max_singleton, cores=cores, plot = plot)
  message("Netboost: Finished clustering step.")
  
  message("Netboost: Finished Netboost.")
  return(results)
}

#' Calculate network adjacencies for filter
#' 
#' @param datan     Dataset
#' @param filter    Filter-Matrix
#' @param softPower Integer. Exponent of the transformation
#' @return Vector with adjacencies for the filter
calculate_adjacency <- function(datan=NULL,filter=NULL,softPower=2) {
  return(sapply(1:nrow(filter),
                function(i) {
                  abs(cor(datan[,filter[i,1]],datan[,filter[i,2]]))^softPower
                }))
}

#' Calculate distance
#' (external wrapper for internal C++ function)
#' Parallelisation inside C++ program with RcppParallel.
#' 
#' @param filter Filter-Matrix
#' @param datan     Dataset
#' @param softPower Integer. Exponent of the transformation
#' @param cores Integer. Amount of CPU cores used (<=1 : sequential)
#' @param verbose   Logical. Additional diagnostic messages.
#' @return Vector with distances (same length as adjacency)
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
#' @param filter    Filter-Matrix
#' @param dist      Filter-Matrix
#' @param max_singleton     The maximal singleton in the clustering. Usually equals the number of features.
#' @param verbose   Logical. Additional diagnostic messages.
#' @param cores     Integer. Amount of CPU cores used.
#' @return Dendrogram
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
#' @param forest Matrix
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
#' @param datan Dataset
#' @param forest Matrix
#' @return List of tree specific objects including dendrogram, tree data and features.
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
  
#' Module detection for an individual tree
#'
#' @param tree_dendro List of tree specific objects including dendrogram, tree data and features originating from the tree_dendro function.
#' @param datan Dataset
#' @param MEDissThres Module Eigengene Dissimilarity Threshold for merging close modules.
#' @param plot Logical. Create plot
#' @param minClusterSize XXX
#' @param name_of_tree XXX
#' @return Object of class hclust
cut_dendro <- function(tree_dendro=NULL, minClusterSize= 10L, 
                       datan=NULL, MEDissThres = NULL,
                       name_of_tree="", plot = TRUE) {
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
  return(list(colors=mergedColors,MEs=MEs,varExplained=MEList$varExplained))
}

#' Module detection for the results from a nb_mcupgma call
#'
#' @param trees List of trees, where one tree is a list of ids and rows
#' @param plot Logical. Create plot.
#' @param datan XXX
#' @param forest XXX
#' @param minClusterSize XXX
#' @param MEDissThres XXX
#' @return List
cut_trees <- function(trees=NULL, datan=NULL, 
                      forest=NULL, minClusterSize= 10L,
                      MEDissThres = NULL, plot = TRUE) {
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
    res[[i]][["varExplained"]] <- cut_dendro$varExplained
    i <- i+1
  }
  return(res)
}

#' Netboost clustering step
#'
#' @param filter    Filter-Matrix as returned by nb_filter.
#' @param dist      Distance matrix as returned by nb_dist.
#' @param datan Dataset
#' @param max_singleton   Integer. The maximal singleton in the clustering. Usually equals the number of features.
#' @param minClusterSize  Integer. The minimum number of features in one module.
#' @param MEDissThres Numeric. Module Eigengene Dissimilarity Threshold for merging close modules.
#' @param cores Integer. Amount of CPU cores used (<=1 : sequential)
#' @param plot Logical. Create plot.
#' @return List
#' @export
nb_clust <- function(filter = NULL,
                     dist = NULL,
                     datan = NULL,
                     max_singleton = dim(datan)[2],
                     minClusterSize = 10L,
                     MEDissThres = 0.25,
                     cores = getOption("mc.cores", 2L),
                     plot = TRUE) {
  forest <- nb_mcupgma(filter=filter,dist=dist,max_singleton=max_singleton,cores=cores)
  trees <- tree_search(forest)
  results <- cut_trees(trees=trees,datan=datan, forest=forest, minClusterSize = minClusterSize, MEDissThres = MEDissThres, plot = plot)
  sum_res <- nb_summary(clust_res = results, plot = plot)
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
        colnames(tmp_MEs_new)[colnames(tmp_MEs)==paste0("ME",col)] <- paste0("ME", (n_MEs))
      }
      if (col == 0 & length(unique(tmp.col)) != 1) {
        n_MEs_background <- n_MEs_background + 1
        tmp.col.new[tmp.col == col] <- -n_MEs_background
        colnames(tmp_MEs_new)[colnames(tmp_MEs)=="ME0"] <- paste0("ME0_", n_MEs_background)
      }
    }
    res$colors <- c(res$colors, tmp.col.new)
    if("MEs" %in% names(res)){res$MEs <- cbind(res$MEs, tmp_MEs_new)}else{res$MEs <- tmp_MEs_new}
    if("varExplained" %in% names(res)){res$varExplained <- c(res$varExplained , clust_res[[tree]]$varExplained)}else{res$varExplained <-  clust_res[[tree]]$varExplained}
  }
  names(res$varExplained)<- colnames(res$MEs)
  
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
    plot_colors[plot_colors <= 0] <- 0
    for (tree in 1:length(res$dendros)) {
      par(mar = c(0, 4, 8, 4))
      plot(res$dendro[[tree]],labels=FALSE)
      par(mar = c(4, 4, 0, 4))
      first_col <- last_col + 1
      last_col <- last_col + length(res$dendro[[tree]]$labels)

      
      plotColorUnderTree(res$dendro[[tree]], c(gray(level = 0.7),
                                               colorspace::rainbow_hcl(n = (length(
                                                 unique(plot_colors)
                                               ) - 1)))[plot_colors[first_col:last_col] + 1])
      # rainbow_hcl(n = (length(unique(plot_colors))-1)))[plot_colors[first_col:last_col]+1])
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
  
  if(!identical(sort(nb_summary$names), sort(colnames(new_data)))){
    stop("The features in new_data (colnames) must correspond to the features in nb_summary (nb_summary$names).")
  }
  
  new_data <- new_data[,nb_summary$names]
  MEs <- moduleEigengenes(new_data, colors = nb_summary$colors)$eigengenes

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

#' TCGA covariates for 180 AML patients.
#'
#' @format A data frame with 188 rows and 78 variables:
#' @source \url{http://www.tcga.com/}
"tcga_aml_covariates"