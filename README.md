# Netboost
Boosting supported network analysis for high-dimensional omics applications.

This package comes bundled with the MC-UPGMA clustering package by Yaniv Loewenstein.

# Requirements
This package is only working on MacOS and Linux (no Windows atm).

Required for building are C/C++ compilers, GNU make, GZIP, Perl.

# Installation
```
# Skip if already installed
#install.packages("devtools")      # Install from Github

# Install netboost
devtools::install_github("PascalSchlosser/netboost")
```

# Example
```R
vignette("netboost")
```
