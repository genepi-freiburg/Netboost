# Netboost
Boosting supported network analysis for high-dimensional omics applications.

This package comes bundled with the MC-UPGMA clustering package by Yaniv Loewenstein.

# Requirements
This package is only working on MacOS and Linux (no Windows atm).

Required for building are C/C++ compilers, GNU make, GZIP, Perl.

# Installation
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("netboost", version = "3.9")
```

# Example
```R
browseVignettes("netboost")
```
