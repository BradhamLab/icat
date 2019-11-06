#! /bin/bash -l
#$ -N seurat_test
#$ -M dyh0110@bu.edu
#$ -m eas

export PYTHONPATH="{$PYTHONPATH}:/projectnb/bradham/PythonModules"
export R_LIBS="{$R_LIBS}:/projectnb/bradham/RPackages"
source activate icat
snakemake --cluster \
'qsub -v PYTHONPATH={cluster.python} -v LD_LIBRARY_PATH={cluster.gcc} -P {cluster.project} -pe omp 4 -e {cluster.error} -o {cluster.out} -l cpu_arch={cluster.cpus}' --jobs 100 --latency-wait 30 --cluster-config cluster.json
