#! /bin/bash -l
#$ -N scanorama_test
#$ -M dyh0110@bu.edu
#$ -m eas

export PYTHONPATH="{$PYTHONPATH}:/projectnb/bradham/PythonModules"
export R_LIBS="{$R_LIBS}:/projectnb/bradham/RPackages"
source activate icat
snakemake --cluster \
'qsub -v PYTHONPATH=/projectnb/bradham/PythonModules -v LD_LIBRARY_PATH=/share/pkg.7/gcc/8.1.0/install/lib64:/share/pkg.7/gcc/8.1.0/install/lib -P bradham -pe omp 4 -e snake_error -o snake_out' --jobs 100 --latency-wait 30 --cluster-config cluster-config.json
