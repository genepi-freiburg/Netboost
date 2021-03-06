% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netboost.R
\name{cut_dendro}
\alias{cut_dendro}
\title{Module detection for an individual tree}
\usage{
cut_dendro(tree_dendro, min_cluster_size = 2L, datan, ME_diss_thres,
  name_of_tree = "", qc_plot = TRUE, n_pc = 1,
  nb_min_varExpl = 0.5)
}
\arguments{
\item{tree_dendro}{List of tree specific objects including dendrogram, 
tree
data and features originating from the tree_dendro function.}

\item{min_cluster_size}{Integer. The minimum number of features in one 
module.}

\item{datan}{Data frame were rows correspond to samples and columns to
features.}

\item{ME_diss_thres}{Numeric. Module Eigengene Dissimilarity Threshold for
merging close modules.}

\item{name_of_tree}{String. Annotating plots and messages.}

\item{qc_plot}{Logical. Should plots be created?}

\item{n_pc}{Number of principal components and variance explained
entries to be calculated. The number of returned variance explained entries
is currently ‘min(n_pc,10)’. If given ‘n_pc’ is greater than 10, a warning 
is issued.}

\item{nb_min_varExpl}{Minimum proportion of variance explained for
returned module eigengenes. The number of PCs is capped at n_pc.}
}
\value{
List
}
\description{
Module detection for an individual tree
}
