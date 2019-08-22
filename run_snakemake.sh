#! /bin/bash -l
#$ -N evaluate_icat
#$ -M dyh0110@bu.edu
#$ -m eas

export PYTHONPATH="{PYTHONPATH}:/projectnb/bradham/PythonModules"
export R_LIBS="{R_LIBS}:/projectnb/bradham/RPackages"
source activate icat
snakemake --cluster 'qsub -v PYTHONPATH=/projectnb/bradham/PythonModules -P bradham -pe omp 3' --jobs 100
