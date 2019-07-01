#! /bin/bash -l
#$ -N evaluate_icat
#$ -M dyh0110@bu.edu
#$ -m eas

export PYTHONPATH="{PTYHONPATH}:/projectnb/bradham/PythonModules"
source activate icat
snakemake --cluster 'qsub -v PYTHONPATH=/projectnb/bradham/PythonModules -P bradham -pe omp 3' --jobs 100
