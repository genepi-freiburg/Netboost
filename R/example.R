# load data
# methylation and RNA data
# load(file.path(netboostPackagePath(), "data/CHR18.Rdata"))
data(tcga_aml_meth_rna_chr18)
data(tcga_aml_covariates)


# filter
nb_filter_boosting(datan=tcga_aml_meth_rna_chr18, stepno=20L, until=0L, progress=1000L, cores=1L, mode=2L)

# adjacencies

# dist

# mcupgma

# h clust

