conda create -n "liana" python=3.9 ipython
conda activate "liana"
conda install -y -c conda-forge mamba
mamba install -y -c conda-forge r-rjson=0.2.21 r-devtools=2.4.3 r-cairo=1.5-14 r-biocmanager=1.30.16
mamba install -y -c bioconda bioconductor-biobase=2.54.0 bioconductor-biocneighbors=1.12.0
Rscript setup_liana_env.r
# mamba install -y -c anaconda jupyter nb_conda_kernels
# mamba install -c conda-forge r-irkernel
