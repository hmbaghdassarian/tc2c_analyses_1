#!/bin/bash
# take in environment name: taken from https://www.shellscript.sh/tips/getopt/
NAME="cellchat"
usage()
{
  echo "Usage: bash -i setup_cellchat_env.sh [ -n | --name ENV_NAME ]"
  exit 2
}

PARSED_ARGUMENTS=$(getopt -o n: --long name: -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
fi

eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -n | --name) NAME="$2" ; shift 2 ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1 - this should not happen."
       usage ;;
  esac
done

# initialize and check conda environment
# echo "NAME   : $NAME"
conda create -y -n "$NAME"
conda activate "$NAME"
conda info|egrep "conda version|active env"

ACT_ENV="$(conda info|egrep "active environment")"
ACT_ENV=(${ACT_ENV// : / })
ACT_ENV=${ACT_ENV[2]}  
if [[ "$ACT_ENV" != "$NAME" ]]; then
  echo "The environment $NAME has not been activated"
  usage
fi

# begin package installs
conda install -y -c conda-forge mamba=0.12.2
mamba install -y -c bioconda bioconductor-biobase=2.50.0
mamba install -y -c conda-forge r-devtools=2.4.0
mamba install -y -c conda-forge r-cairo=1.5_12.2 r-biocmanager=1.30.12  
# mama install -y -c conda-forge r-hdf5r # for SeuratDisk only
Rscript setup_cellchat_env.r
mamba install -y -c conda-forge r-rhpcblasctl=0.20_137
mamba install -y -c conda-forge r-docopt=0.7.1
mamba install -y -c r r-stringr=1.4.0
echo "Complete, activate environment using: conda activate $NAME"
