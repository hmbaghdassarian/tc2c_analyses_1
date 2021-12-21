library(devtools)
devtools::install_github("renozao/repotools@0f76e52253c08063084074061f8ca2c05e8a4818", quiet = T, upgrade_dependencies = F)
# source('http://renozao.github.io/repotools/install.R') #1.10.1

library(repotools)
repotools::install.pkgs('NMF', force = F, quiet = T) #>=0.23.0, should be automatic from repotools

devtools::install_github("jokergoo/circlize@feccef4f5f25f0ad389ccb6d3e36a2e6813ad7a8", quiet = T, upgrade_dependencies = F)

library(BiocManager)
BiocManager::install("ComplexHeatmap", ask = F, update = F) #2.6.2 - should go with appropriate bioconductor version

devtools::install_github("sqjin/CellChat@b3ccf96664e29e702e0c23d0bb4e4fc566c72034")

BiocManager::install("rhdf5", ask = F, update = F) # version 2.34.0

#library(remotes)
#remotes::install_github("mojaveazure/seurat-disk@163f1aade5bac38ed1e9e9c912283a7e74781610") # for SeuratDisk only