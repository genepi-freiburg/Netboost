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

# Contact
If you have any issues using the package then please get in touch with Pascal Schlosser (pascal.schlosser at uniklinik-freiburg.de).
Bug reports etc are most welcome, we want the package to be easy to use for everyone!