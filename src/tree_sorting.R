library(compiler)
enableJIT(3)

library(parallel)
source("/home/jo/Quellcodes/R.tools/tools.R")

## Execute R program as comparison?
COMPARE_WITH_R <- FALSE
BUILD_DENDROGRAM <- FALSE

tic()   # Time-measurement overall.

load("data.Rdata")

tic("load tree")
netboost_forest <- read.table(file = "./DistTOM.mcupgma_tree",
                              row.names = NULL,
                              col.names = c("cluster_id1", "cluster_id2",
                                            "distance", "cluster_id3"),
                              nrow = -1)   ## 100000 for testing
netboost_forest <- as.matrix(netboost_forest)
toc("load tree")

##****************************************************************************
## Replace jo tree-sort
##****************************************************************************
library(Rcpp)

## C++ program setup.
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")
sourceCpp(file.path(".", "rcpp_tree_sort.cpp"))

tic("c++")

## Produce exactly the same list than "trees" in the original (below)
## R program.
## (one entry per tree, each tree is again a list with two vectors
## "ids" and "rows").
## The C++ program delivers both "ids" and "rows" unique and sorted.
rescpp <- rcpp_tree_search(netboost_forest);

toc("c++", echo=TRUE)

##stop("C++ only, no comparison with R.")
##****************************************************************************
##****************************************************************************

if (COMPARE_WITH_R) {
    ## Execute R program as comparison (for step by step).

    ### sort into the disjunct trees ###
    trees <- vector("list", 1)
    trees[[1]]$ids <- as.vector(netboost_forest[1, c(1, 2, 4)])
    trees[[1]]$rows <- c(1)
    sorted.ids <- trees[[1]]$ids

    s <- proc.time()[3]

    tic("loop")
    for (i in 2:nrow(netboost_forest)) {
        if (i %% 1000 == 0) {
            cat("Step: ", i, "/", format(digits = 2,
                                         i * 100 / nrow(netboost_forest)), "% ",
                (proc.time()[3] - s), "\n"
                )
        }

        ## Current row
        new.ids <- as.vector(netboost_forest[i, c(1, 2, 4)])

        ## Kamen die IDs schon mal irgendwo vor?
        hits <- sum(new.ids %in% sorted.ids)

        if (hits == 0) {
            tic("hit 0")
            trees[[length(trees) + 1]] <- list(ids = c(), rows = c())
            trees[[length(trees)]]$ids <- new.ids
            trees[[length(trees)]]$rows <- c(i)
            toc("hit 0")
            #        cat("R: NEW TREE: ", length(trees), "\n")
        } else if (hits == 1) {
            tic("hit 1")
            for (j in 1:length(trees)){
                if (any(new.ids %in% trees[[j]]$ids)){
                    trees[[j]]$ids <- unique(c(trees[[j]]$ids, new.ids))
                    trees[[j]]$rows <- c(trees[[j]]$rows, i)
                    #                cat("R: ADD TO TREE: ", j, "\n")
                    break
                }
            }
            toc("hit 1")
        } else {   ## CHG: 2/3 Hits zusammengefasst
            tic("hit 2/3")
            merge <- c()
            for (j in 1:length(trees)) {
                if (any(new.ids %in% trees[[j]]$ids)) {
                    merge <- c(merge, j)
                    if (length(merge) == hits) {
                        break
                    }
                }
            }
            toc("hit 2/3")

            tic("merge in loop")

            ## Alle Teilbäume in merge[1] zusammenführen.
            for (k in rev(merge[-1])) {
                trees[[merge[1]]]$ids <- c(trees[[merge[1]]]$ids, trees[[k]]$ids)
                trees[[merge[1]]]$rows <- c(trees[[merge[1]]]$rows, trees[[k]]$rows)
                #            cat("R: MERGE TREES: ", k, " < ", merge[[1]], "\n")
                trees[[k]] <- NULL
            }

            trees[[merge[1]]]$ids <- unique(c(trees[[merge[1]]]$ids, new.ids))
            trees[[merge[1]]]$rows <- c(trees[[merge[1]]]$rows, i)

            toc("merge in loop")
        }

        sorted.ids <- c(sorted.ids, new.ids)
        ##    print(sorted.ids)
        #    print(paste("R: ROW: ", i, " TREES: ", length(trees)))
    }
    toc("loop")

    ## Save intermediate result
    save(file="RESULT_trees.RData", trees)

    ##****************************************************************************
    ## COMPARISON C++ <-> R results.
    ##****************************************************************************
    
    res <- rep(FALSE, length(trees))

    print(paste("R Trees  :", length(trees)))
    print(paste("C++ Trees:", length(rescpp)))

    for (i in 1:length(trees)) {
        trees[[i]]$ids <- unique(sort(trees[[i]]$ids));
        trees[[i]]$rows <- unique(sort(trees[[i]]$rows));
    }
    for (i in 1:length(rescpp)) {
        rescpp[[i]]$ids <- unique(sort(rescpp[[i]]$ids));
        rescpp[[i]]$rows <- unique(sort(rescpp[[i]]$rows));
    }

    cpp_result <- rescpp

    for (i in 1:length(trees)) {
        ## Check if row is in rescpp result (list may have different order...)
        for (x in 1:length(rescpp)) {
            if (isTRUE(all.equal(trees[[i]]$ids, rescpp[[x]]$ids)) &&
                isTRUE(all.equal(trees[[i]]$rows, rescpp[[x]]$rows))) {
                res[i] <- TRUE
                
                rescpp[[x]] <- NULL
                
                break
            }
        }

        if (!res[i]) {
            cat("No match on row: ", i, " (R list)\n");
        }
    }

    cat("IDENTICAL ROWS: ", sum(res), " SHOULD BE: ", length(trees), "\n")

    if (sum(res) == length(trees)) {
        print("n'kay!")
    } else {
        print("PROBLEM!")
    }

##    stop("End comparison")
}

print("finished. Now: Dendrogram")

##****************************************************************************
## From here, original program.
## Further input (trees) coming from the C++ program now (rescpp)
##****************************************************************************
trees <- rescpp


tic("postproc")
if (length(trees) > 1){
    ### select the largest tree for analysis ###
    index <- which.max(sapply(sapply(trees, "[[", 1), length))

    index.features <- trees[[index]]$ids[trees[[index]]$ids <= dim(data)[2]]

    colnames.tree <- colnames(data)[index.features]
    cat("WARNING: The largest tree of a forest was selected.",(dim(data)[2]-length(colnames.tree)),"features were removed.\n")
    netboost_tree <- netboost_forest[trees[[index]]$rows,]
    #adjust for removed features
    merge <- sapply(trees, "[[", 1)
    merge[[index]] <- NULL
    merge <- rev(sort(unique(unlist(merge))))
    for (i in merge){
        for (j in c(1, 2, 4)){
            netboost_tree[netboost_tree[,j]>i, j] <- netboost_tree[netboost_tree[,j] > i, j] - 1
        }
    }
    rm.data <- data[, -index.features]
    data <- data[, index.features]
    index.rm.features <- (1:dim(data)[2])[-index.features]
    feature.names.netboost <- colnames(data)
    ## rm(index,index.names)
    rm(index)
} else {
    ### select whole tree ###
    feature.names.netboost <- colnames(data)
    netboost_tree <- netboost_forest
    colnames.tree <- colnames(data)
}

### Dendrogram ###
if (BUILD_DENDROGRAM) {
    cutpoint <- dim(data)[2]
    dendro.netboost <- list()
    dendro.netboost$merge <- netboost_tree[, c(1, 2)]
    dendro.netboost$merge[dendro.netboost$merge <= cutpoint] <- -dendro.netboost$merge[dendro.netboost$merge <= cutpoint]
    dendro.netboost$merge[dendro.netboost$merge > 0] <- dendro.netboost$merge[dendro.netboost$merge > 0] - cutpoint
    dendro.netboost$merge <-  apply(dendro.netboost$merge, c(1, 2), function(x) (as.integer(x))) # why do I need to transform?
    dendro.netboost$height <- netboost_tree[, 3]
    dendro.netboost$order <- 1:cutpoint
    dendro.netboost$labels <- colnames.tree
    class(dendro.netboost) <- "hclust"
    b <- as.dendrogram(dendro.netboost)
    b.order <- order.dendrogram(b)
    dendro.netboost$order <- b.order

    save(file="RESULT_dendronetboost.Rdata", dendro.netboost)

    #pdf(file="dendrogram_netboost.pdf", height = 10, width = 500)
    plot(dendro.netboost, xlab = '', sub = '', lab = FALSE, hang = 0.01)
    dev.off()
}
toc("postproc")

toc()

tictoc_sum("default")
