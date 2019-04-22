## Formatting changed by RStudio-Autoformatting (*mostly* compatible to BioC
## preferences)

## Load WGCNA, try to hide the welcome-message (only a try as it is printed...)
## Workaround using environment to force WGCNA skipping it's welcome-message.
Sys.setenv(ALLOW_WGCNA_THREADS = 1)
suppressPackageStartupMessages(require(WGCNA))
Sys.unsetenv("ALLOW_WGCNA_THREADS")

## require(colorspace)
## require(parallel)

#' Netboost clustering.
#'
#' The Netboost clustering is performed in three subsequent steps. First, a
#' filter of important edges in the network is calculated. Next, pairwise
#' distances are calculated. Last, clustering is performed. For details see
#' Schlosser et al. doi...
#'
#' @name netboost
#' @param datan     Data frame were rows correspond to samples and columns to
#'   features.
#' @param stepno    Integer amount of boosting steps applied in the filtering
#'   step
#' @param until     Stop at index/column (if 0: iterate through all columns).
#'   For testing purposes in large datasets.
#' @param progress  Integer. If > 0, print progress after every X steps (
#' Progress might not be reported 100% accurate due to parallel execution)
#' @param mode      Integer. Mode (0: x86, 1: FMA, 2: AVX). Features are only
#'   available if compiled accordingly and available on the hardware.
#' @param soft_power Integer. Exponent of the transformation. Set automatically
#'   based on the scale free topology criterion if unspecified.
#' @param max_singleton   Integer. The maximal singleton in the clustering.
#'   Usually equals the number of features.
#' @param plot      Logical. Should plots be created?
#' @param min_cluster_size  Integer. The minimum number of features in one module.
#' @param ME_diss_thres Numeric. Module Eigengene Dissimilarity Threshold for
#'   merging close modules.
#' @param cores     Integer. Amount of CPU cores used (<=1 : sequential)
#' @param scale     Logical. Should data be scaled and centered?
#' @param method	A character string specifying the method to be used for
#'   correlation coefficients.
#' @param verbose   Additional diagnostic messages.
#' @param n_pc        Number of principal components and variance explained
#'   entries to be calculated. The number of returned variance explained entries
#'   is currently ‘min(n_pc,10)’. If given ‘n_pc’ is greater than 10, a warning is
#'   issued.
#' @param nb_min_varExpl        Minimum proportion of variance explained for
#'   returned module eigengenes. The number of PCs is capped at n_pc.
#' @return dendros  A list of dendrograms. For each fully separate part of the
#'   network an individual dendrogram.
#' @return names    A vector of feature names.
#' @return colors   A vector of numeric color coding in matching order of names
#'   and module eigengene names (color = 3 -> variable in ME3).
#' @return MEs      Aggregated module measures (Module eigengenes).
#' @return var_explained  Proportion of variance explained per module eigengene
#'   per principal component (max n_pc principal components are listed).
#' @return rotation Matrix of variable loadings divided by their singular
#'   values. datan %*% rotation = MEs (with datan potentially scaled)
#'
#'
#' @examples
#' data('tcga_aml_meth_rna_chr18',  package='netboost')
#' results <- netboost(datan=tcga_aml_meth_rna_chr18, stepno=20L,
#'    soft_power=3L, min_cluster_size=10L, n_pc=2, scale=TRUE,
#'    ME_diss_thres=0.25, plot=TRUE)
#'
#' @export
netboost <-
    function(datan = NULL,
             stepno = 20L,
             until = 0L,
             progress = 1000L,
             mode = 2L,
             soft_power = NULL,
             max_singleton = ncol(datan),
             plot = TRUE,
             min_cluster_size = 2L,
             ME_diss_thres = 0.25,
             n_pc = 1,
             nb_min_varExpl = 0.5,
             cores = as.integer(getOption("mc.cores", 2)),
             scale = TRUE,
             method = c("pearson", "kendall", "spearman"),
             verbose = getOption("verbose")) {
        # Initialize parallelization of WGCNA package.
        if (cores > 1)
            WGCNA::allowWGCNAThreads(nThreads = as.numeric(cores))
        
        if (is.null(datan) || !is.data.frame(datan))
            stop("netboost: Error: datan must be a data.frame.")

        if (is.null(stepno) || !is.integer(stepno))
            stop("netboost: Error: stepno must be an integer.")

        if (is.null(min_cluster_size) || !is.integer(min_cluster_size))
            stop("netboost: Error: min_cluster_size must be an integer.")

        if (is.null(ME_diss_thres) || !is.numeric(ME_diss_thres))
            stop("netboost: Error: ME_diss_thres must be numeric.")

        if (ncol(datan) > 5e+06) {
            stop(
                    "A bug in sparse UPGMA currently prevents analyses",
                    " with more than 5 million features."
            )
        }
        
        if (scale) {
            if(verbose>=1){message("Netboost: Scaling and centering data.")}
            datan <-
                as.data.frame(scale(datan, center = TRUE, scale = TRUE))
        }
        
        if(verbose>=1){message("Netboost: Initialising filter step.")}
        filter <-
            nb_filter(
                datan = datan,
                stepno = stepno,
                until = until,
                progress = progress,
                cores = cores,
                mode = mode
            )
        
        if(verbose>=1){message("Netboost: Finished filter step.")}
        
        if (is.null(soft_power)) {
            # Random subset out of allocation
            random_features <-
                sample(ncol(datan), min(c(10000, ncol(datan))))
            # Call the network topology analysis function
            sft <- WGCNA::pickSoftThreshold(datan[, random_features])
            soft_power <- sft[["powerEstimate"]]
            if(verbose>=0){message(
                paste0(
                    "Netboost: soft_power was set to ",
                    soft_power,
                    " based on the scale free topology criterion."
                )
            )}
        }
        
        if(verbose>=1){message("Netboost: Initialising distance calculation.")}
        dist <-
            nb_dist(
                datan = datan,
                filter = filter,
                soft_power = soft_power,
                cores = cores,
                method = method
            )
        if(verbose>=1){message("Netboost: Finished distance calculation.")}
        
        if(verbose>=1){message("Netboost: Initialising clustering step.")}
        results <-
            nb_clust(
                datan = datan,
                filter = filter,
                dist = dist,
                min_cluster_size = min_cluster_size,
                ME_diss_thres = ME_diss_thres,
                n_pc = n_pc,
                nb_min_varExpl = nb_min_varExpl,
                max_singleton = max_singleton,
                cores = cores,
                plot = plot
            )
        if(verbose>=1){message("Netboost: Finished clustering step.")}
        
        if(verbose>=1){message("Netboost: Finished Netboost.")}
        invisible(results)
    }

#' Calculate network adjacencies for filter
#'
#' @name calculate_adjacency
#' @param datan     Data frame were rows correspond to samples and columns to
#'   features.
#' @param filter    Filter-Matrix as generated by the nb_filter function.
#' @param soft_power Integer. Exponent of the transformation. Set automatically
#'   based on the scale free topology criterion if unspecified.
#' @param method	A character string specifying the method to be used for
#'   correlation coefficients.
#' @return Vector with adjacencies for the filter
calculate_adjacency <-
    function(datan,
             filter,
             soft_power = 2,
             method = c("pearson", "kendall", "spearman")) {
        return(vapply(X=seq(
            from = 1,
            to = nrow(filter),
            by = 1
        ),
        FUN=function(i) {
            abs(WGCNA::cor(datan[, filter[i, 1]],
                           datan[, filter[i, 2]], method = method)) ^ soft_power
        },
        FUN.VALUE=1))
    }

#' Calculate distance (external wrapper for internal C++ function)
#' Parallelisation inside C++ program with RcppParallel.
#'
#' @name nb_dist
#' @param filter    Filter-Matrix as generated by the nb_filter function.
#' @param datan     Data frame were rows correspond to samples and columns to
#'   features.
#' @param soft_power Integer. Exponent of the transformation. Set automatically
#'   based on the scale free topology criterion if unspecified.
#' @param cores     Integer. Amount of CPU cores used (<=1 : sequential).
#' @param verbose   Additional diagnostic messages.
#' @param method	A character string specifying the method to be used for
#'   correlation coefficients.
#' @return Vector with distances.
#'
#' @examples
#'  data('tcga_aml_meth_rna_chr18', package='netboost')
#'  cores <- as.integer(getOption('mc.cores', 2))
#'  datan <- as.data.frame(scale(tcga_aml_meth_rna_chr18, center=TRUE,
#'  scale=TRUE))
#'  filter <- nb_filter(datan=datan, stepno=20L, until=0L, progress=1000L,
#'  cores=cores,mode=2L)
#'  dist <- nb_dist(datan=datan, filter=filter, soft_power=3L, cores=cores)
#'  summary(dist)
#'
#' @export
nb_dist <-
    function(filter,
             datan,
             soft_power = 2,
             cores = getOption("mc.cores", 2L),
             verbose = getOption("verbose"),
             method = c("pearson", "kendall", "spearman")) {
        # if (is.null(filter) || is.null(adjacency)) stop('Both filter and
        # adjacency must be provided')
        
        if (!(is.matrix(filter) &&
              (nrow(filter) > 0) && (ncol(filter) > 0)))
            stop("filter must be matrix with nrow > 0 and ncol >0")
        
        # if (!(is.vector(adjacency) && (length(adjacency) > 0)))
        # stop('adjacency is required a vector with length > 0')
        cores <- max(cores, 1)
        
        ## RcppParallel amount of threads to be started
        setThreadOptions(numThreads = cores)
        
        return(cpp_dist_tom(
            filter,
            calculate_adjacency(
                datan = datan,
                filter = filter,
                soft_power = soft_power,
                method = method)
        ))
    }

#' Calculate dendrogram for a sparse distance matrix
#' (external wrapper MC-UPGMA clustering package Loewenstein et al.
#'
#' @name nb_mcupgma
#' @param filter Filter-Matrix as generated by the nb_filter function.
#' @param dist Distance-Matrix as generated by the nb_dist function.
#' @param max_singleton Integer The maximal singleton in the clustering.
#'                      Usually equals the number of features.
#' @param cores Integer Amount of CPU cores used (<=1 : sequential)
#' @param verbose Logical Additional diagnostic messages.
#' @return Raw dendrogram to be processed by tree_search and tree_dendro.
#'
#' @examples
#'    data('tcga_aml_meth_rna_chr18', package='netboost')
#'    cores <- as.integer(getOption('mc.cores', 2))
#'    datan <- as.data.frame(scale(tcga_aml_meth_rna_chr18,
#'    center=TRUE, scale=TRUE))
#'    filter <- nb_filter(datan=datan, stepno=20L, until=0L,
#'                        progress=1000L, cores=cores, mode=2L)
#'    dist <- nb_dist(datan=datan, filter=filter, soft_power=3L, cores=cores)
#'    max_singleton = dim(tcga_aml_meth_rna_chr18)[2]
#'    forest <- nb_mcupgma(filter=filter, dist=dist,
#'                         max_singleton=max_singleton, cores=cores)
#'  head(forest)
#'
#' @export
nb_mcupgma <-
    function(filter,
             dist,
             max_singleton,
             cores = getOption("mc.cores",
                               2L),
             verbose = getOption("verbose")) {
        # Deletes all files under netboostTmpPath(), esp. clustering/iteration_
        netboostTmpCleanup()
        
        if (max_singleton > 5e+06) {
            stop("A bug in sparse UPGMA currently prevents analyses",
                       " with more than 5 million features.")
        }
        
        if (!dir.create(file.path(netboostTmpPath(), "clustering")))
            stop(
                "Unable to create: ",
                file.path(netboostTmpPath(), "clustering")
            )
        
        file_dist_edges <-
            file.path(netboostTmpPath(), "clustering", "dist.edges")
        file_dist_tree <-
            file.path(netboostTmpPath(), "clustering", "dist.mcupgma_tree")
        
        # write.table(file='clustering/dist.edges',
        write.table(
            file = file_dist_edges,
            cbind(
                format(filter, scientific = FALSE,
                       trim = TRUE),
                format(dist, scientific = FALSE, trim = TRUE)
            ),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
        )
        
        file_dist_edges_gz <- paste(file_dist_edges, ".gz", sep="")
        ret <- R.utils::gzip(file_dist_edges, destname = file_dist_edges_gz,
                             overwrite = TRUE)

        if (attr(ret, "nbrOfBytes") <= 0 || !R.utils::isGzipped(file_dist_edges_gz))
            warning(paste("Gzip maybe failed on:", file_dist_edges,
                          "Return:", as.character(ret),
                          " Bytes: ", attr(ret, "nbrOfBytes")))

        ret <-
            mcupgma_exec(
                exec = "cluster.pl",
                "-max_distance",
                1,
                "-max_singleton",
                max_singleton,
                "-iterations 1000 -heap_size 10000000 -num_hash_buckets 40",
                "-jobs",
                cores,
                "-retries 1",
                "-output_tree_file",
                file_dist_tree,
                "-split_unmodified_edges",
                max(cores,2L),
                file_dist_edges_gz,
                console = FALSE
            )
        
        if (verbose)
            message(ret)
        
        if (!file.exists(file_dist_tree) ||
            file.info(file_dist_tree)[["size"]] == 0)
            stop("No output file created. mcupgma error.")
        
        return(as.matrix(
            read.table(
                file = file_dist_tree,
                row.names = NULL,
                col.names = c("cluster_id1",
                              "cluster_id2", "distance", "cluster_id3")
            )
        ))
    }

#' Extracts independent trees from nb_mcupgma results (external wrapper for
#' internal C++ function)
#'
#' @name tree_search
#' @param forest Raw dendrogram-matrix as generated by the nb_mcupgma function.
#' @return List
#'   
#' @examples
#'    data('tcga_aml_meth_rna_chr18',  package='netboost')
#'    cores <- as.integer(getOption('mc.cores', 2))
#'    datan <- as.data.frame(scale(tcga_aml_meth_rna_chr18, center=TRUE,
#'                                 scale=TRUE))
#'    filter <- nb_filter(datan=datan, stepno=20L, until=0L, progress=1000L,
#'                        cores=cores,mode=2L)
#'    dist <- nb_dist(datan=datan, filter=filter, soft_power=3L, cores=cores)
#'    max_singleton = dim(tcga_aml_meth_rna_chr18)[2]
#'    forest <- nb_mcupgma(filter=filter, dist=dist,
#'                         max_singleton=max_singleton, cores=cores)
#'    trees <- tree_search(forest)
#'    str(trees[[length(trees)]])
#'
#' @export
tree_search <- function(forest = NULL) {
    ## Check for integer values cannot be done here, as either the user must
    ## have set up all values with as.integer() or R delivers default numeric
    ## (double). In that case, Rcpp converts the matrix.  (Alternative: convert
    ## manually 'matrix(as.integer(forest), nrow=nrow(forest))')
    #forest <- as.matrix(forest)
    
    if (is.null(forest) || !is.matrix(forest))
        stop("forest must be provided (as matrix)")
    
    return(cpp_tree_search(forest))
}


#' Calculate the dendrogram for an individual tree
#'
#' @name tree_dendro
#' @param tree A list with two elements. ids, which is an integer vector of
#'   feature identifiers and rows, which is an integer vector of selected rows
#'   in the corresponding forest
#' @param datan     Data frame were rows correspond to samples and columns to
#'   features.
#' @param forest Raw dendrogram-matrix as generated by the nb_mcupgma function.
#' @return List of tree specific objects including dendrogram, tree data and
#'   features.
tree_dendro <- function(tree,
                        datan,
                        forest) {
    index.features <- tree[["ids"]][tree[["ids"]] <= dim(datan)[2]]
    data_tree <- datan[, index.features]
    
    colnames_tree <- colnames(datan)[index.features]
    tree_cluster <- forest[tree[["rows"]], , drop = FALSE]
    
    all_ids <- rev(sort(unique(c(forest[, c(1, 2, 4)]))))
    none_tree_ids <- all_ids[!(all_ids %in% tree[["ids"]])]
    for (i in none_tree_ids) {
        for (j in c(1, 2, 4)) {
            tree_cluster[tree_cluster[, j] > i, j] <-
                tree_cluster[tree_cluster[,
                                          j] > i, j] - 1
        }
    }
    
    feature_names <- colnames(data_tree)
    
    cutpoint <- dim(data_tree)[2]
    dendro <- list()
    dendro[["merge"]] <- tree_cluster[, c(1, 2), drop = FALSE]
    dendro[["merge"]][dendro[["merge"]] <= cutpoint] <-
        -dendro[["merge"]][dendro[["merge"]] <= cutpoint]
    dendro[["merge"]][dendro[["merge"]] > 0] <-
        dendro[["merge"]][dendro[["merge"]] > 0] - cutpoint
    dendro[["merge"]] <- apply(dendro[["merge"]], c(1, 2), function(x) {
        (as.integer(x))
    })
    dendro[["height"]] <- tree_cluster[, 3]
    dendro[["order"]] <- seq(from = 1, to = cutpoint, by = 1)
    dendro[["labels"]] <- colnames_tree
    class(dendro) <- "hclust"
    b <- as.dendrogram(dendro)
    b.order <- order.dendrogram(b)
    dendro[["order"]] <- b.order
    
    return(list(
        dendro = dendro,
        data = data_tree,
        names = feature_names
    ))
}

#' Module detection for an individual tree
#'
#' @name cut_dendro
#' @param tree_dendro List of tree specific objects including dendrogram, tree
#'   data and features originating from the tree_dendro function.
#' @param min_cluster_size  Integer. The minimum number of features in one module.
#' @param datan     Data frame were rows correspond to samples and columns to
#'   features.
#' @param ME_diss_thres Numeric. Module Eigengene Dissimilarity Threshold for
#'   merging close modules.
#' @param name_of_tree String. Annotating plots and messages.
#' @param plot      Logical. Should plots be created?
#' @param n_pc        Number of principal components and variance explained
#'   entries to be calculated. The number of returned variance explained entries
#'   is currently ‘min(n_pc,10)’. If given ‘n_pc’ is greater than 10, a warning is
#'   issued.
#' @param nb_min_varExpl        Minimum proportion of variance explained for
#'   returned module eigengenes. The number of PCs is capped at n_pc.
#' @return List
cut_dendro <-
    function(tree_dendro,
             min_cluster_size = 2L,
             datan,
             ME_diss_thres,
             name_of_tree = "",
             plot = TRUE,
             n_pc = 1,
             nb_min_varExpl = 0.5) {
        dynamicMods <-
            dynamicTreeCut::cutreeDynamic(
                dendro = tree_dendro[["dendro"]],
                method = "tree",
                deepSplit = TRUE,
                minClusterSize = min_cluster_size
            )
        ### Merging of Dynamic Modules ### Calculate eigengenes
        MEList <-
            netboost::nb_moduleEigengenes(
                expr = tree_dendro[["data"]],
                colors = dynamicMods,
                n_pc = n_pc,
                nb_min_varExpl = nb_min_varExpl
            )
        MEs <- MEList[["nb_eigengenes"]]
        # Calculate dissimilarity of module eigengenes
        MEDiss <- 1 - cor(MEs)
        # Cluster module eigengenes
        if (length(MEDiss) > 1) {
            METree <- hclust(as.dist(MEDiss), method = "average")
            if (plot == TRUE) {
                graphics::plot(
                    METree,
                    main = paste0(name_of_tree,
                                  "Clustering of module eigengenes"),
                    xlab = "",
                    sub = ""
                )
                graphics::abline(h = ME_diss_thres, col = "red")
            }
            
            merged <-
                WGCNA::mergeCloseModules(
                    exprData = tree_dendro[["data"]],
                    dynamicMods,
                    cutHeight = ME_diss_thres,
                    verbose = 3
                )
            mergedColors <- merged[["colors"]]
            # Calculate eigengenes
            MEList <-
                netboost::nb_moduleEigengenes(
                    expr = tree_dendro[["data"]],
                    colors = merged[["colors"]],
                    n_pc = n_pc,
                    nb_min_varExpl = nb_min_varExpl
                )
            MEs <- MEList[["nb_eigengenes"]]
            MEDiss <- 1 - cor(MEs)
            if (length(MEDiss) > 1) {
                METree <- hclust(as.dist(MEDiss), method = "average")
                if (plot == TRUE &
                    length(tree_dendro[["dendro"]][["labels"]]) > 2) {
                    graphics::plot(
                        METree,
                        main = paste0(
                            name_of_tree,
                            "Clustering of merged module eigengenes"
                        ),
                        xlab = "",
                        sub = ""
                    )
                    WGCNA::plotDendroAndColors(
                        dendro = tree_dendro[["dendro"]],
                        colors = mergedColors,
                        "Merged Dynamic",
                        dendroLabels = FALSE,
                        hang = 0.01,
                        addGuide = TRUE,
                        guideHang = 0.05,
                        main = paste0(name_of_tree, "Cluster Dendrogram")
                    )
                }
            }
        } else {
            message("\nOnly one module in ", name_of_tree, ".\n")
            mergedColors <- dynamicMods
            if (length(tree_dendro[["dendro"]][["labels"]]) > 2) {
                if (plot == TRUE) {
                    graphics::plot(
                        tree_dendro[["dendro"]],
                        main = paste0(
                            name_of_tree,
                            "Cluster Dendrogram ",
                            "(Tree maximaly consists out of one module.)"
                        )
                    )
                }
            } else {
                message(
                    "\nOnly two elements in the one module in ",
                    name_of_tree,
                    " (no plot generated).\n"
                )
            }
        }
        message(
            "\nNetboost extracted ",
            length(table(mergedColors)),
            " modules (including background) with an average size of ",
            mean(table(mergedColors)[-1]),
            " (excluding background) from ",
            substr(name_of_tree,
                   start = 1, stop = (nchar(name_of_tree) - 1)),
            ".\n"
        )
        return(
            list(
                colors = mergedColors,
                MEs = MEs,
                var_explained = MEList[["var_explained"]],
                rotation = MEList[["rotation"]]
            )
        )
        # svd_PCs = MEList[["eigengenes"]],
    }

#' Module detection for the results from a nb_mcupgma call
#'
#' @name cut_trees
#' @param trees List of trees, where one tree is a list of ids and rows
#' @param datan     Data frame were rows correspond to samples and columns to
#'   features.
#' @param forest Raw dendrogram-matrix as generated by the nb_mcupgma function.
#' @param min_cluster_size  Integer. The minimum number of features in one module.
#' @param ME_diss_thres Numeric. Module Eigengene Dissimilarity Threshold for
#'   merging close modules.
#' @param plot      Logical. Should plots be created?
#' @param n_pc        Number of principal components and variance explained
#'   entries to be calculated. The number of returned variance explained entries
#'   is currently ‘min(n_pc,10)’. If given ‘n_pc’ is greater than 10, a warning is
#'   issued.
#' @param nb_min_varExpl        Minimum proportion of variance explained for
#'   returned module eigengenes. The number of PCs is capped at n_pc.
#' @return List
#'
#' @examples
#' data('tcga_aml_meth_rna_chr18',  package='netboost')
#'  cores <- as.integer(getOption('mc.cores', 2))
#'  datan <- as.data.frame(scale(tcga_aml_meth_rna_chr18, center=TRUE,
#'  scale=TRUE))
#'  filter <- nb_filter(datan=datan, stepno=20L, until=0L, progress=1000L,
#'  cores=cores,mode=2L)
#'  dist <- nb_dist(datan=datan, filter=filter, soft_power=3L, cores=cores)
#'  max_singleton = dim(tcga_aml_meth_rna_chr18)[2]
#'  forest <- nb_mcupgma(filter=filter, dist=dist, max_singleton=max_singleton,
#'  cores=cores)
#'  trees <- tree_search(forest)
#'  results <- cut_trees(trees=trees,datan=datan, forest=forest,
#'  min_cluster_size=10L, ME_diss_thres=0.25, plot=TRUE)
#'
#' @export
cut_trees <-
    function(trees,
             datan,
             forest,
             min_cluster_size = 2L,
             ME_diss_thres,
             plot = TRUE,
             n_pc = 1,
             nb_min_varExpl = 0.5) {
        res <- list()
        i <- 1L
        
        for (tree in trees) {
            tree_dendro_res <-
                tree_dendro(tree = tree,
                            datan = datan,
                            forest = forest)
            res[[i]] <- list()
            res[[i]][["dendro"]] <- tree_dendro_res[["dendro"]]
            res[[i]][["data"]] <- tree_dendro_res[["data"]]
            res[[i]][["names"]] <- tree_dendro_res[["names"]]
            cut_dendro_res <-
                cut_dendro(
                    tree_dendro = tree_dendro_res,
                    min_cluster_size = min_cluster_size,
                    datan = datan,
                    ME_diss_thres = ME_diss_thres,
                    name_of_tree = paste0("Tree ",
                                          i, ":"),
                    plot = plot,
                    n_pc = n_pc,
                    nb_min_varExpl = nb_min_varExpl
                )
            res[[i]][["colors"]] <- cut_dendro_res[["colors"]]
            res[[i]][["MEs"]] <- cut_dendro_res[["MEs"]]
            # res[[i]][['svd_PCs']] <- cut_dendro_res[["svd_PCs"]]
            res[[i]][["var_explained"]] <-
                cut_dendro_res[["var_explained"]]
            res[[i]][["rotation"]] <- cut_dendro_res[["rotation"]]
            i <- i + 1
        }
        
        return(res)
    }

#' Netboost clustering step
#'
#' @name nb_clust
#' @param filter    Filter-Matrix as generated by the nb_filter function.
#' @param dist      Distance-Matrix as generated by the nb_dist function.
#' @param datan     Data frame were rows correspond to samples and columns to
#'   features.
#' @param max_singleton   Integer. The maximal singleton in the clustering.
#'   Usually equals the number of features.
#' @param min_cluster_size  Integer. The minimum number of features in one module.
#' @param ME_diss_thres Numeric. Module Eigengene Dissimilarity Threshold for
#'   merging close modules.
#' @param cores     Integer. Amount of CPU cores used (<=1 : sequential)
#' @param plot Logical. Create plot.
#' @param n_pc        Number of principal components and variance explained
#'   entries to be calculated. The number of returned variance explained entries
#'   is currently ‘min(n_pc,10)’. If given ‘n_pc’ is greater than 10, a warning is
#'   issued.
#' @param nb_min_varExpl        Minimum proportion of variance explained for
#'   returned module eigengenes. The number of PCs is capped at n_pc.
#' @return List
#'
#' @examples
#'  data('tcga_aml_meth_rna_chr18',  package='netboost')
#'  cores <- as.integer(getOption('mc.cores', 2))
#'  datan <- as.data.frame(scale(tcga_aml_meth_rna_chr18, center=TRUE, 
#'  scale=TRUE))
#'  filter <- nb_filter(datan=datan, stepno=20L, until=0L, progress=1000L,
#'  cores=cores,mode=2L)
#'  dist <- nb_dist(datan=datan, filter=filter, soft_power=3L, cores=cores)
#'  max_singleton = dim(tcga_aml_meth_rna_chr18)[2]
#'  pdf("test.pdf",width=30)
#'  sum_res <- nb_clust(filter=filter, dist=dist, datan=datan,
#'  max_singleton=max_singleton, min_cluster_size=10L, ME_diss_thres=0.25,
#'  cores=cores, plot=TRUE, n_pc=2L, nb_min_varExpl=0.5)
#'  dev.off()
#' @export
nb_clust <-
    function(filter,
             dist,
             datan,
             max_singleton = dim(datan)[2],
             min_cluster_size = 2L,
             ME_diss_thres = 0.25,
             cores = getOption("mc.cores", 2L),
             plot = TRUE,
             n_pc = 1,
             nb_min_varExpl = 0.5) {
        forest <-
            nb_mcupgma(
                filter = filter,
                dist = dist,
                max_singleton = max_singleton,
                cores = cores
            )
        trees <- tree_search(forest)
        results <-
            cut_trees(
                trees = trees,
                datan = datan,
                forest = forest,
                min_cluster_size = min_cluster_size,
                ME_diss_thres = ME_diss_thres,
                plot = plot,
                n_pc = n_pc,
                nb_min_varExpl = nb_min_varExpl
            )
        sum_res <- nb_summary(clust_res = results, plot = plot)
        
        return(sum_res)
    }

#' Summarize results from a forest. Plot trees together.
#'
#' @name nb_summary
#' @param clust_res Clustering results from cut_trees call.
#' @param plot Logical. Create plot.
#' @return List
#'
#' @examples
#' data('tcga_aml_meth_rna_chr18',  package='netboost')
#'  cores <- as.integer(getOption('mc.cores', 2))
#'  datan <- as.data.frame(scale(tcga_aml_meth_rna_chr18, center=TRUE,
#'  scale=TRUE))
#'  filter <- nb_filter(datan=datan, stepno=20L, until=0L, progress=1000L,
#'  cores=cores,mode=2L)
#'  dist <- nb_dist(datan=datan, filter=filter, soft_power=3L, cores=cores)
#'  max_singleton = dim(tcga_aml_meth_rna_chr18)[2]
#'  forest <- nb_mcupgma(filter=filter,dist=dist,max_singleton=max_singleton,
#'  cores=cores)
#'  trees <- tree_search(forest)
#'  results <- cut_trees(trees=trees,datan=datan, forest=forest,
#'  min_cluster_size=10L, ME_diss_thres=0.25, plot=FALSE)
#'  sum_res <- nb_summary(clust_res=results, plot=TRUE)
#'
#' @export
nb_summary <- function(clust_res, plot = TRUE) {
    res <- vector("list")
    # res$clust_res <- clust_res
    n_MEs <- 0
    n_MEs_background <- 0
    for (tree in seq(from = 1,
                     to = length(clust_res),
                     by = 1)) {
        res[["dendros"]][[tree]] <- clust_res[[tree]][["dendro"]]
        res[["names"]] <- c(res[["names"]], clust_res[[tree]][["names"]])
        tmp.col <- clust_res[[tree]][["colors"]]
        tmp.col.new <- tmp.col
        tmp_MEs <- clust_res[[tree]][["MEs"]]
        tmp_MEs_new <- tmp_MEs
        tmp_rotation <- clust_res[[tree]][["rotation"]]
        tmp_rotation <-
            do.call("cbind", lapply(
                tmp_rotation,
                FUN = function(x) {
                    y <-
                        matrix(0,
                               nrow = length(clust_res[[tree]][["names"]]),
                               ncol = ncol(x))
                    rownames(y) <- clust_res[[tree]][["names"]]
                    y[rownames(x),] <- x
                    colnames(y) <- colnames(x)
                    return(y)
                }
            ))
        tmp_rotation_new <- tmp_rotation
        for (col in unique(tmp.col)) {
            if (col != 0 | length(unique(tmp.col)) == 1) {
                n_MEs <- n_MEs + 1
                tmp.col.new[tmp.col == col] <- n_MEs
                colnames(tmp_MEs_new)[grepl(pattern = paste0("ME", col, "_"),
                                            colnames(tmp_MEs))] <-
                    gsub(
                        pattern = paste0("ME",
                                         col, "_"),
                        replacement = paste0("ME", (n_MEs), "_"),
                        colnames(tmp_MEs_new)[grepl(pattern = paste0("ME",
                                                                     col, "_"),
                                                    colnames(tmp_MEs))]
                    )
                colnames(tmp_rotation_new)[grepl(pattern = paste0("ME", col,
                                                                  "_"),
                                                 colnames(tmp_rotation))] <-
                    gsub(
                        pattern = paste0("ME", col,
                                         "_"),
                        replacement = paste0("ME", (n_MEs), "_"),
                        colnames(tmp_rotation_new)[grepl(pattern =
                                                paste0("ME", col, "_"),
                                                colnames(tmp_rotation))]
                    )
            }
            if (col == 0 & length(unique(tmp.col)) != 1) {
                n_MEs_background <- n_MEs_background + 1
                tmp.col.new[tmp.col == col] <- -n_MEs_background
                colnames(tmp_MEs_new)[grepl(pattern = "ME0_",
                                            colnames(tmp_MEs))] <-
                    gsub(
                        pattern = "ME0_",
                        replacement = paste0("ME0_", n_MEs_background, "_"),
                        colnames(tmp_MEs_new)[grepl(pattern = "ME0_",
                                                    colnames(tmp_MEs))]
                    )
                colnames(tmp_rotation_new)[grepl(pattern = "ME0_",
                                                 colnames(tmp_rotation))] <-
                    gsub(
                        pattern = "ME0_",
                        replacement = paste0("ME0_", n_MEs_background, "_"),
                        colnames(tmp_rotation_new)[grepl(pattern = "ME0_",
                                                     colnames(tmp_rotation))]
                    )
            }
        }
        res[["colors"]] <- c(res[["colors"]], tmp.col.new)
        if ("MEs" %in% names(res)) {
            res[["MEs"]] <- cbind(res[["MEs"]], tmp_MEs_new)
        } else {
            res[["MEs"]] <- tmp_MEs_new
        }
        if ("rotation" %in% names(res)) {
            res[["rotation"]] <- c(res[["rotation"]], list(tmp_rotation_new))
        } else {
            res[["rotation"]] <- list(tmp_rotation_new)
        }
        ## if('svd_PCs' %in% names(res)){res$svd_PCs <- cbind(res$svd_PCs,
        ## clust_res[[tree]]$svd_PCs)}else{res$svd_PCs <-
        ## clust_res[[tree]]$svd_PCs}
        if ("var_explained" %in% names(res)) {
            res[["var_explained"]] <-
                cbind(res[["var_explained"]], clust_res[[tree]][["var_explained"]])
        } else {
            res[["var_explained"]] <- clust_res[[tree]][["var_explained"]]
        }
    }
    rownames(res[["var_explained"]]) <-
        paste0("PC", seq(
            from = 1,
            to = nrow(res[["var_explained"]]),
            by = 1
        ))
    colnames(res[["var_explained"]]) <-
        unique(unlist(lapply(
            strsplit(split = "_pc",
                     colnames(res[["MEs"]])),
            FUN = function(x) {
                x[1]
            }
        )))
    
    res[["rotation"]] <-
        do.call("cbind", lapply(
            res[["rotation"]],
            FUN = function(x) {
                y <- matrix(0, nrow = length(res[["names"]]), ncol = ncol(x))
                rownames(y) <- res[["names"]]
                y[rownames(x),] <- x
                colnames(y) <- colnames(x)
                return(y)
            }
        ))
    
    message(
        "\nNetboost detected ",
        n_MEs,
        " modules and ",
        n_MEs_background,
        " background modules in ",
        length(clust_res),
        " trees resulting in ",
        ncol(res[["MEs"]]),
        " aggreagate measures.\n"
    )
    message("Average size of the modules was ",
        mean(table(res[["colors"]][!(res[["colors"]] <= 0)])), ".\n")
    message(
        sum(res[["colors"]] <= 0),
        " of ",
        length(res[["colors"]]),
        " features (",
        (sum(res[["colors"]] <=
                 0) * 100 /
             length(res[["colors"]])),
        "%) were not assigned to modules.\n"
    )
    
    if (plot == TRUE) {
        nb_plot_dendro(nb_summary = res, labels = FALSE)
    }
    return(res)
}

#' Transfer of Netboost clustering to new data.
#'
#' @name nb_tranfer
#' @param nb_summary Netboost results as generated by the nb_summary function.
#' @param new_data Data frame were rows correspond to samples and columns to
#'   features.
#' @param scale     Logical. Should data be scaled and centered?
#' @param only_module_membership     Logical. Should only module memberships be
#'   transfered and PCs be newly computed?
#' @return List
#'
#' @examples
#' data('tcga_aml_meth_rna_chr18',  package='netboost')
#' results <- netboost(datan = tcga_aml_meth_rna_chr18, stepno = 20L,
#'     soft_power = 3L, min_cluster_size = 10L, n_pc = 2, scale=TRUE,
#'     ME_diss_thres = 0.25, plot=FALSE)
#' ME_transfer <- nb_transfer(nb_summary = results,
#'     new_data = tcga_aml_meth_rna_chr18,
#'     scale = TRUE)
#' all(round(results[["MEs"]], 12) == round(ME_transfer, 12))
#'
#' @export
nb_transfer <-
    function(nb_summary = NULL,
             new_data = NULL,
             scale = FALSE,
             only_module_membership = FALSE) {
        if (!exists("new_data"))
            stop("datan must be provided")
        
        if (!(is.data.frame(new_data) &&
              (nrow(new_data) > 0) && (ncol(new_data) >
                                       0)))
            stop("new_data must be a data frame with dim() > (0,0).")
        
        if (length(nb_summary[["colors"]]) != ncol(new_data)) {
            stop("The number of features in new_data must ",
                       "correspond to the number in nb_summary.")
        }
        
        if (!identical(sort(nb_summary[["names"]]), sort(colnames(new_data)))) {
            stop("The features in new_data (colnames) must ", 
                       "correspond to the features in nb_summary ",
                       "(nb_summary$names).")
        }
        
        new_data <- new_data[, nb_summary[["names"]]]
        
        if (scale) {
            new_data <-
                as.data.frame(scale(new_data, center = TRUE, scale = TRUE))
        }
        
        if (!only_module_membership) {
            MEs <- as.matrix(new_data) %*% nb_summary[["rotation"]]
        } else {
            MEs <- netboost::nb_moduleEigengenes(expr = new_data,
                                                 colors = nb_summary[["colors"]])[["nb_eigengenes"]]
            colnames(MEs)[lapply(strsplit(x = colnames(MEs), split = "-"), FUN = length) >
                              1] <-
                paste0("ME0_", substring(text = colnames(MEs)[lapply(strsplit(x = colnames(MEs),
                                                                              split = "-"), FUN = length) > 1], first = 4))
        }
        
        MEs <- MEs[, colnames(nb_summary[["MEs"]])]
        rownames(MEs) <- rownames(new_data)
        return(MEs)
    }


#' Boosting via C++ function. Parallelisation by R-package parallel with forking
#' (overhead of this method does not fall into account as single steps are
#' ~10s).
#'
#' Parallelisation via multicore (via 'parallel'-package). So *nix only atm.
#'
#' @name nb_filter
#' @param datan     Data frame were rows correspond to samples and columns to
#'   features.
#' @param stepno    Integer amount of boosting steps
#' @param until     Stop at index/column (if 0: iterate through all columns)
#' @param progress  Integer. If > 0, print progress after every X steps (mind:
#'   parallel!)
#' @param cores     Integer. Amount of CPU cores used (<=1 : sequential)
#' @param mode      Integer. Mode (0: x86, 1: FMA, 2: AVX). Features are only
#'   available if compiled accordingly and available on the hardware.
#' @return matrix n times 2 matrix with the indicies of the n unique entrees of
#'   the filter
#'
#'
#' @examples
#' data('tcga_aml_meth_rna_chr18',  package='netboost')
#'  cores <- as.integer(getOption('mc.cores', 2))
#'  datan <- as.data.frame(scale(tcga_aml_meth_rna_chr18, center=TRUE,
#'  scale=TRUE))
#'  filter <- nb_filter(datan=datan, stepno=20L, until=0L, progress=1000L,
#'  cores=cores,mode=2L)
#'  head(filter)
#'  nrow(filter)/(ncol(datan)*(ncol(datan)-1)/2) # proportion of potential undirected edges
#'
#' @export
nb_filter <-
    function(datan,
             stepno = 20L,
             until = 0L,
             progress = 1000L,
             cores = getOption("mc.cores",
                               2L),
             mode = 2L) {
        if (!exists("datan"))
            stop("datan must be provided")
        
        if (!(is.data.frame(datan) &&
              (nrow(datan) > 0) && (ncol(datan) > 0)))
            stop("datan must be a data frame with dim() > (0,0).")
        
        if (!(is.integer(stepno) && (stepno > 0)))
            stop("stepno must be an integer > 0.")
        
        if (!(is.integer(until) && (until >= 0)))
            stop("until must be an integer >= 0.")
        
        if (until == 0)
            until <- ncol(datan)
        
        # check mode
        if (!(mode %in% c(0, 1, 2))) {
            stop("mode must be 0 (x86), 1 (FMA) or 2 (AVX).")
        }
        
        if (ncol(datan) > 5e+06) {
            stop("A bug in sparse UPGMA currently prevents analyses ",
                       "with more than 5 million features.")
        }
        
        message(paste("Netboost: Filtering"))
        
        ## Initialize data structures for optimized boosting (once)
        cpp_filter_base(as.matrix(datan), stepno, mode_ = mode)
        
        ## Parallelization 'conventional' via mclapply.
        if (cores > 1) {
            # print(paste('Parallel version:', cores, 'cores'))
            
            boosting_filter <- mclapply(seq(1, until), function(x) {
                if ((((x - 1) %% progress) == 0)) {
                    message(sprintf("idx: %d (%.1f%%) - %s", x,
                                  x * 100 / until, date()))
                }
                
                cpp_filter_step(x)
            }, mc.cores = cores)
        } else {
            ## Sequential function for debugging.  print(paste('Sequential version'))
            boosting_filter <- lapply(seq(1, until), function(x) {
                if ((((x - 1) %% progress) == 0)) {
                    message(sprintf("idx: %d (%.1f%%) - %s", x,
                                  x * 100 / until, date()))
                }
                
                cpp_filter_step(x)
            })
        }
        
        ## Important!: stop (free memory, else suitable memory is still blocked)
        cpp_filter_end()
        
        filter <-
            do.call("rbind", lapply(seq_along(boosting_filter), function(x) {
                return(as.data.frame(cbind(
                    as.integer(boosting_filter[[x]]),
                    as.integer(rep(x,
                                   length(
                                       boosting_filter[[x]]
                                   )))
                )))
            }))

        filter <- unique(t(apply(filter, 1, sort)))
        colnames(filter) <- c("cluster_id1", "cluster_id2")
        rownames(filter) <- seq(from = 1,
                                to = nrow(filter),
                                by = 1)

        return(filter)
    }

#' Plot dendrogram from Netboost output.
#'
#' @name nb_plot_dendro
#' @param nb_summary Netboost results as generated by the nb_summary function.
#' @param labels TRUE/FALSE indicator of whether labels should be attached to
#'   the leafs.
#' @param main Plot title.
#' @param colorsrandom TRUE/FALSE indicator of whether module colors should be
#'   shuffeled.
#' @return invisible null
#'
#' @examples
#' data('tcga_aml_meth_rna_chr18',  package='netboost')
#' results <- netboost(datan = tcga_aml_meth_rna_chr18, stepno = 20L,
#' soft_power = 3L, min_cluster_size = 10L, n_pc = 2, scale=TRUE,
#' ME_diss_thres = 0.25, plot = FALSE)
#' set.seed(1234) # reproducible but shuffled color-module matching
#' nb_plot_dendro(nb_summary = results, labels = FALSE, main = 'Test',
#' colorsrandom = TRUE)
#'
#' @export
nb_plot_dendro <-
    function(nb_summary = NULL,
             labels = FALSE,
             main = "",
             colorsrandom = FALSE) {
        if (!exists("nb_summary"))
            stop("Netboost output (nb_summary) must be provided.")
        
        colorHeight <- 0.2
        graphics::layout(matrix(seq(
            from = 1,
            to = (2 * length(nb_summary[["dendros"]])),
            by = 1
        ),
        nrow = 2),
        heights = c(1 - colorHeight, colorHeight))
        
        last_col <- 0
        n_colors <- length(unique(nb_summary[["colors"]]))
        middle <- floor(n_colors / 2)
        if (colorsrandom) {
            shuffel_index <- sample(x = n_colors, size = n_colors)
        } else {
            shuffel_index <- seq(from = 1,
                                 to = n_colors,
                                 by = 1)
        }
        shuffel_index <-
            c(shuffel_index[(2 * middle + 1)][(2 * middle + 1) == n_colors],
              rbind(shuffel_index[seq(from = 1, to = middle, by = 1)], shuffel_index[seq(from = (middle +
                                                                                                     1),
                                                                                         to = (2 * middle),
                                                                                         by = 1)]))
        plot_colors <-
            colorspace::rainbow_hcl(n = (length(unique(
                nb_summary[["colors"]]
            ))))[shuffel_index][as.factor(nb_summary[["colors"]])]
        plot_colors[nb_summary[["colors"]] <= 0] <- grDevices::gray(level = 0.7)
        for (tree in seq_along(nb_summary[["dendros"]])) {
            graphics::par(mar = c(0, 4, 8, 4))
            first_col <- last_col + 1
            last_col <-
                last_col + length(nb_summary[["dendros"]][[tree]][["labels"]])
            if (labels) {
                graphics::plot(nb_summary[["dendros"]][[tree]],
                     labels = nb_summary[["names"]][seq(from = first_col,
                                                   to = last_col,
                                                   by = 1)],
                     main = main)
            } else {
                graphics::plot(nb_summary[["dendros"]][[tree]],
                     labels = FALSE,
                     main = main)
            }
            graphics::par(mar = c(4, 4, 0, 4))
            WGCNA::plotColorUnderTree(nb_summary[["dendros"]][[tree]],
                                      colors = plot_colors[seq(from = first_col,
                                                               to = last_col,
                                                               by = 1)],
                                      rowLabels = "")
        }
        invisible(NULL)     # Only to not have empty value section
    }


#' Netboost module aggregate extraction.
#'
#' This is a modification of WGCNA::moduleEigengenes() (version WGCNA_1.66) to
#' include more than the first principal component. For details see
#' WGCNA::moduleEigengenes().
#'
#' @name nb_moduleEigengenes
#' @param expr     Expression data for a single set in the form of a data frame
#'   where rows are samples and columns are genes (probes).
#' @param colors    A vector of the same length as the number of probes in
#'   ‘expr’, giving module color for all probes (genes). Color ‘'grey'’ is
#'   reserved for unassigned genes.     Expression
#' @param n_pc       Number of principal components and variance explained
#'   entries to be calculated. The number of returned variance explained entries
#'   is currently ‘min(n_pc,10)’. If given ‘n_pc’ is greater than 10, a warning is
#'   issued.
#' @param align     Controls whether eigengenes, whose orientation is
#'   undetermined, should be aligned with average expression (‘align = 'along
#'   average'’, the default) or left as they are (‘align = ''’). Any other value
#'   will trigger an error.
#' @param exclude_grey   Should the improper module consisting of 'grey' genes be
#'   excluded from the eigengenes?
#' @param grey          Value of ‘colors’ designating the improper module. Note
#'   that if ‘colors’ is a factor of numbers, the default value will be
#'   incorrect.
#' @param subHubs       Controls whether hub genes should be substituted for
#'   missing eigengenes. If ‘TRUE’, each missing eigengene (i.e., eigengene
#'   whose calculation failed and the error was trapped) will be replaced by a
#'   weighted average of the most connected hub genes in the corresponding
#'   module. If this calculation fails, or if ‘subHubs==FALSE’, the value of
#'   ‘trapErrors’ will determine whether the offending module will be removed or
#'   whether the function will issue an error and stop.
#' @param trapErrors    Controls handling of errors from that may arise when
#'   there are too many ‘NA’ entries in expression data. If ‘TRUE’, errors from
#'   calling these functions will be trapped without abnormal exit.  If ‘FALSE’,
#'   errors will cause the function to stop. Note, however, that ‘subHubs’ takes
#'   precedence in the sense that if ‘subHubs==TRUE’ and ‘trapErrors==FALSE’, an
#'   error will be issued only if both the principal component and the hubgene
#'   calculations have failed.
#' @param return_valid_only   logical; controls whether the returned data frame of
#'   module eigengenes contains columns corresponding only to modules whose
#'   eigengenes or hub genes could be calculated correctly (‘TRUE’), or whether
#'   the data frame should have columns for each of the input color labels
#'   (‘FALSE’).
#' @param soft_power     The power used in soft-thresholding the adjacency
#'   matrix. Only used when the hubgene approximation is necessary because the
#'   principal component calculation failed. It must be non-negative. The
#'   default value should only be changed if there is a clear indication that it
#'   leads to incorrect results.
#' @param scale         logical; can be used to turn off scaling of the
#'   expression data before calculating the singular value decomposition. The
#'   scaling should only be turned off if the data has been scaled previously,
#'   in which case the function can run a bit faster. Note however that the
#'   function first imputes, then scales the expression data in each module. If
#'   the expression contain missing data, scaling outside of the function and
#'   letting the function impute missing data may lead to slightly different
#'   results than if the data is scaled within the function.
#' @param verbose       Controls verbosity of printed progress messages. 0 means
#'   silent, up to (about) 5 the verbosity gradually increases.
#' @param indent        A single non-negative integer controlling indentation of
#'   printed messages. 0 means no indentation, each unit above that adds two
#'   spaces.
#' @param nb_min_varExpl        Minimum proportion of variance explained for
#'   returned module eigengenes. Is capped at n_pc.
#'
#' @return eigengenes   Module eigengenes in a dataframe, with each column
#'   corresponding to one eigengene. The columns are named by the corresponding
#'   color with an ‘'ME'’ prepended, e.g., ‘MEturquoise’ etc. If
#'   ‘return_valid_only==FALSE’, module eigengenes whose calculation failed have
#'   all components set to ‘NA’.
#' @return averageExpr  If ‘align == 'along average'’, a dataframe containing
#'   average normalized expression in each module. The columns are named by the
#'   corresponding color with an ‘'AE'’ prepended, e.g., ‘AEturquoise’ etc.
#' @return var_explained A dataframe in which each column corresponds to a
#'   module, with the component ‘var_explained[PC, module]’ giving the variance
#'   of module ‘module’ explained by the principal component no. ‘PC’. The
#'   calculation is exact irrespective of the number of computed principal
#'   components. At most 10 variance explained values are recorded in this
#'   dataframe.
#' @return n_pc          A copy of the input ‘n_pc’.
#' @return validMEs     A boolean vector. Each component (corresponding to the
#'   columns in ‘data’) is ‘TRUE’ if the corresponding eigengene is valid, and
#'   ‘FALSE’ if it is invalid. Valid eigengenes include both principal
#'   components and their hubgene approximations. When ‘return_valid_only==FALSE’,
#'   by definition all returned eigengenes are valid and the entries of
#'   ‘validMEs’ are all ‘TRUE’.
#' @return validColors  A copy of the input colors with entries corresponding to
#'   invalid modules set to ‘grey’ if given, otherwise 0 if ‘colors’ is numeric
#'   and 'grey' otherwise.
#' @return allOK        Boolean flag signalling whether all eigengenes have been
#'   calculated correctly, either as principal components or as the hubgene
#'   average approximation.
#' @return allPC        Boolean flag signalling whether all returned eigengenes
#'   are principal components.
#' @return isPC         Boolean vector. Each component (corresponding to the
#'   columns in ‘eigengenes’) is ‘TRUE’ if the corresponding eigengene is the
#'   first principal component and ‘FALSE’ if it is the hubgene approximation or
#'   is invalid.
#' @return isHub        Boolean vector. Each component (corresponding to the
#'   columns in ‘eigengenes’) is ‘TRUE’ if the corresponding eigengene is the
#'   hubgene approximation and ‘FALSE’ if it is the first principal component or
#'   is invalid.
#' @return validAEs     Boolean vector. Each component (corresponding to the
#'   columns in ‘eigengenes’) is ‘TRUE’ if the corresponding module average
#'   expression is valid.
#' @return allAEOK      Boolean flag signalling whether all returned module
#'   average expressions contain valid data. Note that ‘return_valid_only==TRUE’
#'   does not imply ‘allAEOK==TRUE’: some invalid average expressions may be
#'   returned if their corresponding eigengenes have been calculated correctly.
#' @export
nb_moduleEigengenes <-
    function(expr,
             colors,
             n_pc = 1,
             align = "along average",
             exclude_grey = FALSE,
             grey = if (is.numeric(colors))
                 0
             else
                 "grey",
             subHubs = TRUE,
             trapErrors = FALSE,
             return_valid_only = trapErrors,
             soft_power = 6,
             scale = TRUE,
             verbose = 0,
             indent = 0,
             nb_min_varExpl = 0.5) {
        spaces <- indentSpaces(indent)
        
        if (verbose == 1)
            message(
                paste(
                    spaces,
                    "moduleEigengenes: Calculating",
                    nlevels(as.factor(colors)),
                    "module eigengenes in given set."
                )
            )
        
        if (is.null(expr)) {
            stop("moduleEigengenes: Error: expr is NULL. ")
        }
        
        if (is.null(colors)) {
            stop("moduleEigengenes: Error: colors is NULL. ")
        }
        
        if (is.null(dim(expr)) || length(dim(expr)) != 2)
            stop("moduleEigengenes: Error: expr must be two-dimensional.")
        
        if (dim(expr)[2] != length(colors))
            stop(
                "moduleEigengenes: Error: ncol(expr) and length(colors) must be equal (one color per gene)."
            )
        
        if (is.factor(colors)) {
            nl <- nlevels(colors)
            nlDrop <- nlevels(colors[, drop = TRUE])
            if (nl > nlDrop)
                stop(
                        "Argument 'colors' contains unused levels (empty modules). ",
                        "Use colors[, drop=TRUE] to get rid of them."
                )
        }
        
        if (soft_power < 0)
            stop("soft_power must be non-negative")
        
        alignRecognizedValues <- c("", "along average")
        
        if (!is.element(align, alignRecognizedValues)) {
            message(
                paste(
                    "ModulePrincipalComponents: Error:",
                    "parameter align has an unrecognised value:",
                    align,
                    "; Recognized values are ",
                    alignRecognizedValues
                )
            )
            stop()
        }
        
        maxVarExplained <- 10
        
        if (n_pc > maxVarExplained)
            warning(paste("Given n_pc is too large. Will use value", maxVarExplained))
        
        nVarExplained <- min(n_pc, maxVarExplained)
        modlevels <- levels(factor(colors))
        
        if (exclude_grey)
            if (sum(as.character(modlevels) != as.character(grey)) > 0) {
                modlevels <-
                    modlevels[as.character(modlevels) != as.character(grey)]
            } else {
                stop(
                        "Color levels are empty. Possible reason: the only color is grey",
                        " and grey module is excluded from the calculation."
                )
            }
        
        PrinComps <- data.frame(matrix(NA, nrow = dim(expr)[[1]],
                                       ncol = length(modlevels)))
        nb_PrinComps <- data.frame(matrix(NA, nrow = dim(expr)[[1]],
                                          ncol = 0))
        rotation <- vector("list", length(modlevels))
        averExpr <-
            data.frame(matrix(NA, nrow = dim(expr)[[1]],
                              ncol = length(modlevels)))
        varExpl <-
            data.frame(matrix(NA, nrow = nVarExplained,
                              ncol = length(modlevels)))
        validMEs <- rep(TRUE, length(modlevels))
        validAEs <- rep(FALSE, length(modlevels))
        isPC <- rep(TRUE, length(modlevels))
        isHub <- rep(FALSE, length(modlevels))
        validColors <- colors
        names(PrinComps) <-
            paste(moduleColor.getMEprefix(), modlevels, sep = "")
        names(averExpr) <- paste("AE", modlevels, sep = "")
        for (i in seq_along(modlevels)) {
            if (verbose > 1)
                message(paste(
                    spaces,
                    "moduleEigengenes : Working on ME for module",
                    modlevels[i]
                ))
            modulename <- modlevels[i]
            restrict1 <-
                as.character(colors) == as.character(modulename)
            if (verbose > 2)
                message(paste(spaces, " ...", sum(restrict1), "features"))
            datModule <- as.matrix(t(expr[, restrict1]))
            n <- dim(datModule)[1]
            p <- dim(datModule)[2]
            pc <- try({
                if (verbose > 5)
                    message(paste(spaces, " ...scaling"))
                if (scale)
                    datModule <- t(scale(t(datModule)))
                if (verbose > 5)
                    message(paste(spaces, " ...calculating SVD"))
                svd1 <-
                    svd(datModule,
                        nu = min(n, p, n_pc),
                        nv = min(n, p, n_pc))
                nb_PCA <-
                    stats::prcomp(
                        x = t(datModule),
                        retx = TRUE,
                        center = FALSE,
                        scale. = FALSE,
                        tol = NULL,
                        rank. = NULL
                    )
                nb_PCA[["x"]] <- t(t(nb_PCA[["x"]]) / svd1[["d"]])
                nb_PCA[["rotation"]] <- t(t(nb_PCA[["rotation"]]) / svd1[["d"]])
                if (verbose > 5)
                    message(paste(spaces, " ...calculating PVE"))
                veMat <-
                    WGCNA::cor(svd1[["v"]][, seq(from = 1,
                                            to = min(n, p, nVarExplained),
                                            by = 1)], t(datModule), use = "p")
                varExpl[seq(
                    from = 1,
                    to = min(n, p, nVarExplained),
                    by = 1
                ), i] <- rowMeans(veMat ^ 2,
                                  na.rm = TRUE)
                svd1[["v"]][, 1]
            }, silent = TRUE)
            
            if (methods::is(pc, "try-error")) {
                if (!trapErrors)
                    stop(pc)
                if (verbose > 0) {
                    message(
                        paste(
                            spaces,
                            " ..ME calculation of module",
                            modulename,
                            "failed with the following error:"
                        )
                    )
                    message(
                        paste(
                            spaces,
                            "     ",
                            pc,
                            spaces,
                            " ..the offending module has been removed."
                        )
                    )
                }
                warning(
                    paste(
                        "Eigengene calculation of module",
                        modulename,
                        "failed with the following error \n     ",
                        pc,
                        "The offending module has been removed.\n"
                    )
                )
                validMEs[i] <- FALSE
                isPC[i] <- FALSE
                isHub[i] <- FALSE
                validColors[restrict1] <- grey
            } else {
                PrinComps[, i] <- pc
                nb_n_pcs <-
                    min(c(which(
                        cumsum(varExpl[seq(
                            from = 1,
                            to = min(n, p,
                                     nVarExplained),
                            by = 1
                        ), i]) > nb_min_varExpl
                    ), nVarExplained))
                nb_PrinComps <-
                    base::cbind(nb_PrinComps, nb_PCA[["x"]][, seq(from = 1,
                                                       to = nb_n_pcs,
                                                       by = 1)])
                colnames(nb_PrinComps)[seq(
                    from = (ncol(nb_PrinComps) - nb_n_pcs + 1),
                    to = ncol(nb_PrinComps),
                    by = 1
                )] <- paste0(
                    moduleColor.getMEprefix(),
                    modlevels[i],
                    "_pc",
                    seq(
                        from = 1,
                        to = nb_n_pcs,
                        by = 1
                    )
                )
                colnames(nb_PCA[["rotation"]])[seq(from = 1,
                                              to = nb_n_pcs,
                                              by = 1)] <- paste0(
                                                  moduleColor.getMEprefix(),
                                                  modlevels[i],
                                                  "_pc",
                                                  seq(
                                                      from = 1,
                                                      to = nb_n_pcs,
                                                      by = 1
                                                  )
                                              )
                rotation[[i]] <-
                    nb_PCA[["rotation"]][, seq(from = 1,
                                          to = nb_n_pcs,
                                          by = 1),
                                    drop = FALSE]
                ae <- try({
                    if (isPC[i])
                        scaledExpr <- scale(t(datModule))
                    averExpr[, i] <-
                        rowMeans(scaledExpr, na.rm = TRUE)
                    if (align == "along average") {
                        if (verbose > 4)
                            message(
                                paste(
                                    spaces,
                                    " .. aligning module eigengene with average expression."
                                )
                            )
                        corAve <-
                            WGCNA::cor(averExpr[, i], PrinComps[, i], use = "p")
                        if (!is.finite(corAve))
                            corAve <- 0
                        if (corAve < 0)
                            PrinComps[, i] <- -PrinComps[, i]
                    }
                    0
                }, silent = TRUE)
                if (methods::is(ae, "try-error")) {
                    if (!trapErrors)
                        stop(ae)
                    if (verbose > 0) {
                        message(
                            paste(
                                spaces,
                                " ..Average expression calculation of module",
                                modulename,
                                "failed with the following error:"
                            )
                        )
                        message(
                            paste(
                                spaces,
                                "     ",
                                ae,
                                spaces,
                                " ..the returned average expression ",
                                "vector will be invalid."
                            )
                        )
                    }
                    warning(
                        paste(
                            "Average expression calculation of module",
                            modulename,
                            "failed with the following error \n     ",
                            ae,
                            "The returned average expression vector will",
                            " be invalid.\n"
                        )
                    )
                }
                validAEs[i] <- !methods::is(ae, "try-error")
            }
        }
        
        allOK <- (sum(!validMEs) == 0)
        if (return_valid_only && sum(!validMEs) > 0) {
            PrinComps <- PrinComps[, validMEs]
            averExpr <- averExpr[, validMEs]
            varExpl <- varExpl[, validMEs]
            validMEs <- rep(TRUE, times = ncol(PrinComps))
            isPC <- isPC[validMEs]
            isHub <- isHub[validMEs]
            validAEs <- validAEs[validMEs]
        }
        
        allPC <- (sum(!isPC) == 0)
        allAEOK <- (sum(!validAEs) == 0)
        list(
            eigengenes = PrinComps,
            averageExpr = averExpr,
            var_explained = varExpl,
            n_pc = n_pc,
            validMEs = validMEs,
            validColors = validColors,
            allOK = allOK,
            allPC = allPC,
            isPC = isPC,
            isHub = isHub,
            validAEs = validAEs,
            allAEOK = allAEOK,
            nb_eigengenes = nb_PrinComps,
            rotation = rotation
        )
    }

## #' Example to get access to the MCUPGMA executables.
## #' 
## #' @examples
## #' mcupgma_example()
## #' @export
## mcupgma_example <- function() {
##     exec <- netboostMCUPGMAPath()
##     files <- Sys.glob(file.path(exec, "*"))
##     paste("Available MCUPGMA executables and scripts under:", exec)
##     print(sapply(files, basename, USE.NAMES = FALSE))
## }

#' Test/example code. Applies netboost to the TCGA-AML CHR18 DNA methylation and
#' gene expression data supplied with the package.
#'
#' @param cores Integer. CPU cores to use.
#' @param keep Logical. Keep mcupgma intermediate files.
#' @return Netboost result
#'
#' @examples
#' nb_example()
#'
#' @export
nb_example <-
    function(cores = getOption("mc.cores", 2L),
             keep = FALSE) {
        # Keep data local.
        exa_env <- new.env()
        
        # load data methylation and RNA data 180 patients x 5283 features
        data("tcga_aml_meth_rna_chr18",
             package = "netboost",
             envir = exa_env)
        
        pdfFile <- file.path(tempdir(), "results_netboost.pdf")
        
        # pdf(file=file.path(getwd(), 'results_netboost.pdf'), width = 30)
        pdf(file = pdfFile, width = 30)
        results <-
            netboost(
                datan = exa_env[["tcga_aml_meth_rna_chr18"]],
                stepno = 20L,
                soft_power = 3L,
                min_cluster_size = 10L,
                n_pc = 2,
                scale = TRUE,
                ME_diss_thres = 0.25
            )
        # set.seed(1234)
        nb_plot_dendro(nb_summary = results,
                       labels = TRUE,
                       colorsrandom = TRUE)
        dev.off()
        
        if (file.exists(pdfFile)) {
            message(paste0("PDF created:", pdfFile))
            
            # If default PDF viewer is assigned, try to show PDF.
            if (!is.null(getOption("pdfviewer"))) {
      #          system2(getOption("pdfviewer"), pdfFile)
            }
        }
        
        ### Transfer results to the same data (bug check)
        ME_transfer <-
            nb_transfer(
                nb_summary = results,
                new_data = exa_env[["tcga_aml_meth_rna_chr18"]],
                scale = TRUE
            )
        
        all(round(results[["MEs"]], 12) == round(ME_transfer, 12))
        
        # Cleanup all produced temporary filed (esp. clustering/iteration_*)
        if (!keep)
            netboostTmpCleanup()
        else
            message(paste("Kept MCUPGMA temporary files in:", netboostTmpPath()))
        
        invisible(results)
    }
