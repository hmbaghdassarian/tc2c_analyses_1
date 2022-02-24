# install liana
library(devtools)
devtools::install_github('saezlab/liana@976659316aa3574993822eac51bd1b2633bdd0de', quiet = T, upgrade_dependencies = F)
# check here

# note, newer version of CellChat and dependencies than in setup_cellchat_env
# first check which of these are installed
devtools::install_github("renozao/repotools@0f76e52253c08063084074061f8ca2c05e8a4818", quiet = T, upgrade_dependencies = F)
library(repotools)
repotools::install.pkgs('NMF', force = F, quiet = T) #>=0.23.0, should be automatic from repotools
devtools::install_github("jokergoo/circlize@9da95467f0b1a79285dd72e59e04bb0803bce1b7", quiet = T, upgrade_dependencies = F)
library(BiocManager)
BiocManager::install("ComplexHeatmap", ask = F, update = F) #>=2.6.2 - should go with appropriate bioconductor version
devtools::install_github("sqjin/CellChat@55629ca8eee84c9dd6903462505c0cd6dcec12d1", quiet = T, upgrade_dependencies = F)

# update biocparallel and all its dependencies as it casues errors
BiocManager::install("BiocParallel", ask = F, update = T) # version 1.28.3