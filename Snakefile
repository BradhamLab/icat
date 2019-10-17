
import glob
import os
import itertools
# from icat.src import snakemake_utils as utils

configfile: './snakemake-config.yaml'

FILES = ['X.csv', 'obs.csv', 'var.csv']
RUNS = ['Experiment1.Perturbation1Sim1Rep1',
        'Experiment1.Perturbation1Sim1Rep2',
        'Experiment1.Perturbation1Sim1Rep3']
# RUNS = utils.get_simulation_ids(config['simulations']['json'],
#                                 config['simulations']['sims'],
#                                 config['simulations']['reps'])
# EXPERIMENTS = utils.get_experiment_ids(RUNS) 
EXPERIMENTS = ['Experiment1Perturbation1']
METHODS = ['icat', 'seurat233', 'seurat311', 'scanorama', 'icat_scan',
           'seurat_icat']
# METHODS = ['icat', 'seurat233', 'seurat311', 'scanorama']

SIMULATED = ["data/processed/simulated/{run}/{out}".\
             format(run=run, out=out)\
             for run, out in itertools.product(RUNS, FILES)]
BENCHMARK = ['cellmix1', 'cellmix2', 'cellmix3', 'cellmix4', 'sc_celseq2']
MIXES = BENCHMARK[:-1]

rule all:
    input:
        ['data/results/simulated/icat/{run}/performance.csv'.format(run=run)\
          for run in RUNS],
        ['data/results/simulated/seurat233/{run}/obs.csv'.format(run=run)\
          for run in RUNS],
        ['data/results/simulated/scanorama/{run}/obs.csv'.format(run=run)\
          for run in RUNS],
        ['data/results/simulated/{method}/{run}/results.csv'.format(
          method=method, run=run) for method, run in\
          itertools.product(METHODS, RUNS)],
        'data/results/simulated/final/results.csv',
        ['reports/figures/simulated/performance/{exp}_metrics.svg'.format(
               exp=exp) for exp in EXPERIMENTS]
        # 'data/processed/Kang/X.csv',
        # 'data/results/benchmark/results.csv'
        # ['data/processed/benchmark/{bench}/X.csv'.format(bench=bench)\
        #     for bench in BENCHMARK]

# ---------------------------- Process Kang Data -------------------------------

rule format_kang_data:
    params:
        datadir='data/raw/Kang/',
        outdir='data/processed/Kang/'
    output:
        'data/processed/Kang/X.csv',
        'data/processed/Kang/obs.csv',
        'data/processed/Kang/var.csv'
    script:
        'src/generate_kang_et_al.py'
    
# ---------------------------- Generate Simulated Data -------------------------
rule simulate_data:
    input:
        json=config['simulations']['json'],
    params:
        sims=config['simulations']['sims'],
        reps=config['simulations']['reps'],
        outdir="data/processed/simulated/"
    output:
        data=SIMULATED,
        csv="data/processed/simulated/simulations.csv"
    script:
        "src/generate_simulated_datasets.py"

# ---------------------- Fit and Analyze Simulated Data ------------------------
rule fit_simulated:
    input:
        X='data/processed/simulated/{run}/X.csv',
        obs='data/processed/simulated/{run}/obs.csv',
        var='data/processed/simulated/{run}/var.csv'
    params:
        treatment='Treatment',
        control='Control',
        label='Population',
        plotdir='reports/figures/simulated/{run}/'
    output:
        json='data/interim/fits/simulated/{run}_fit.json'
    script:
        'src/fit_louvain.py'

rule simulated_icat:
    input:
        X='data/processed/simulated/{run}/X.csv',
        obs='data/processed/simulated/{run}/obs.csv',
        var='data/processed/simulated/{run}/var.csv',
        json='data/interim/fits/simulated/{run}_fit.json',
        ncfs='data/external/simulated_ncfs_params.json'
    params:
        name='{run}',
        treatment='Treatment',
        control='Control',
        plotdir='reports/figures/simulated/{run}/icat/',
        outdir='data/results/simulated/icat/{run}/',
    output:
        csv='data/results/simulated/icat/{run}/performance.csv',
        obs='data/results/simulated/icat/{run}/obs.csv',
        p1='reports/figures/simulated/{run}/icat/umap_Louvain.svg',
        p2='reports/figures/simulated/{run}/icat/umap_NCFS-Louvain.svg',
        p3='reports/figures/simulated/{run}/icat/umap_NCFS-SSLouvain.svg'
    script:
        'src/evaluate_icat.py'

rule simulated_seurat233:
    input:
        X='data/processed/simulated/{run}/X.csv',
        obs='data/processed/simulated/{run}/obs.csv',
        json='data/interim/fits/simulated/{run}_fit.json'
    params:
        treatment='Treatment',
        control='Control',
        label='Population',
        seurat=config['libraries']['seurat2.3.3']
    output:
        X='data/results/simulated/seurat233/{run}/X.csv',
        obs='data/results/simulated/seurat233/{run}/obs.csv'
    script:
        'src/run_seurat2-3-3.R'

rule simulated_seurat311:
    input:
        X='data/processed/simulated/{run}/X.csv',
        obs='data/processed/simulated/{run}/obs.csv',
        json='data/interim/fits/simulated/{run}_fit.json'
    params:
        treatment='Treatment',
        seurat=config['libraries']['seurat3.1.1']
    output:
        X='data/results/simulated/seurat311/{run}/X.csv',
        obs='data/results/simulated/seurat311/{run}/obs.csv',
    script:
        'src/run_seurat3-1-1.R'

rule simulated_seurat_icat:
    input:
        X='data/results/simulated/seurat311/{run}/X.csv',
        obs='data/results/simulated/seurat311/{run}/obs.csv',
        json='data/interim/fits/simulated/{run}_fit.json',
        ncfs='data/external/simulated_ncfs_params.json'
    params:
        treatment='Treatment',
        controls='Control',
        outdir='data/results/simulated/seurat_icat/{run}/',
        cluster='seurat.sslouvain'
    output:
        X='data/results/simulated/seurat_icat/{run}/X.csv',
        obs='data/results/simulated/seurat_icat/{run}/obs.csv',
    script:
        'src/scanorama_icat.py'

rule simulated_scanorama:
    input:
        X='data/processed/simulated/{run}/X.csv',
        obs='data/processed/simulated/{run}/obs.csv',
        var='data/processed/simulated/{run}/var.csv',
        json='data/interim/fits/simulated/{run}_fit.json',
    output:
        X='data/results/simulated/scanorama/{run}/X.csv',
        obs='data/results/simulated/scanorama/{run}/obs.csv',
        var='data/results/simulated/scanorama/{run}/var.csv'
    params:
        treatment='Treatment',
        controls='Control',
        outdir='data/results/simulated/scanorama/{run}/'
    script:
        'src/run_scanorama.py'

rule simulated_scanorama_icat:
    input:
        X='data/results/simulated/scanorama/{run}/X.csv',
        obs='data/results/simulated/scanorama/{run}/obs.csv',
        var='data/results/simulated/scanorama/{run}/var.csv',
        json='data/interim/fits/simulated/{run}_fit.json',
        ncfs='data/external/simulated_ncfs_params.json'
    output:
        X='data/results/simulated/icat_scan/{run}/X.csv',
        obs='data/results/simulated/icat_scan/{run}/obs.csv',
        var='data/results/simulated/icat_scan/{run}/var.csv'
    params:
        outdir='data/results/simulated/icat_scan/{run}/',
        treatment='Treatment',
        controls ='Control',
        cluster='scanorama.sslouvain'
    script:
        'src/scanorama_icat.py'

rule plot_simulated_runs:
    input:
        X='data/results/simulated/{method}/{run}/X.csv',
        obs='data/results/simulated/{method}/{run}/obs.csv',
        fit='data/interim/fits/simulated/{run}_fit.json',
    output:
        treatment='reports/figures/simulated/clusters/{run}/{method}/treatment_umap.svg',
        known='reports/figures/simulated/clusters/{run}/{method}/known_types_umap.svg',
        cluster='reports/figures/simulated/clusters/{run}/{method}/cluster_umap.svg',
        known_bar='reports/figures/simulated/clusters/{run}/{method}/cluster_distribution.svg',
        cluster_bar='reports/figures/simulated/clusters/{run}/{method}/cell_type_distribution.svg'
    params:
        treatment='Treatment',
        label='Population',
        plotdir='reports/figures/simulated/clusters/{run}/{method}/'
    script:
        'src/plot_simulated.py'

rule evaluate_methods_simulated:
    input:
        obs='data/results/simulated/{method}/{run}/obs.csv'
    output:
        csv='data/results/simulated/{method}/{run}/results.csv'
    params:
        run='{run}',
        method='{method}',
        identity='Population'
    script:
        'src/evaluate_clusters.py'

rule combine_evaluations_simulated:
    input:
        csvs=expand('data/results/simulated/{method}/{run}/results.csv',
                    method=METHODS, run=RUNS)
    output:
        csv='data/results/simulated/final/results.csv'
    script:
        'src/concatenate_results.py'

rule summarize_simulated:
    input:
        results='data/results/simulated/final/results.csv'
    output:
        means=['data/results/simulated/final/{exp}/metric_means.csv'.format(
               exp=exp) for exp in EXPERIMENTS],
        devs=['data/results/simulated/final/{exp}/metric_stds.csv'.format(
               exp=exp) for exp in EXPERIMENTS],
        svgs=['reports/figures/simulated/performance/{exp}_metrics.svg'.format(
               exp=exp) for exp in EXPERIMENTS]
    params:
        plotdir='reports/figures/simulated/performance/',
        outdir='data/results/simulated/final/'
    script:
        'src/summarize_simulated.py'

# --------------------------- Process Benchmark Data ---------------------------

# normalize across cells, subset to highly variable genes, and ln(X +1)
# transform cells
rule format_benchmark_data:
    input:
        counts=['data/raw/benchmark/{bench}.count.csv'.format(bench=bench)\
                for bench in BENCHMARK],
        meta=['data/raw/benchmark/{bench}.metadata.csv'.format(bench=bench)\
              for bench in BENCHMARK]
    params:
        outdir='data/processed/benchmark/',
        plotdir='figures/benchmark/'
    output:
        X='data/processed/benchmark/X.csv',
        obs='data/processed/benchmark/obs.csv',
        var='data/processed/benchmark/var.csv',
        gene_svg='figures/benchmark/filter_genes_dispersion.svg'
    script:
        'src/format_benchmark_data.py'

# just combine separate csvs together, b/c Seurat forces it's own preprocessing
rule combine_benchmark_data:
    input:
        counts=['data/raw/benchmark/{bench}.count.csv'.format(bench=bench)\
                for bench in BENCHMARK],
        meta=['data/raw/benchmark/{bench}.metadata.csv'.format(bench=bench)\
              for bench in BENCHMARK]
    params:
        outdir='data/raw/benchmark/combined/',
    output:
        X='data/raw/benchmark/combined/X.csv',
        obs='data/raw/benchmark/combined/obs.csv',
        var='data/raw/benchmark/combined/var.csv',
    script:
        'src/combine_data.py'
    
# ------------------------ Fit and Analyze Benchmark Data ----------------------
rule fit_benchmark_data:
    input:
        X='data/processed/benchmark/X.csv',
        obs='data/processed/benchmark/obs.csv',
        var='data/processed/benchmark/var.csv' 
    params:
        treatment='benchmark',
        control='sc_celseq2',
        label='mixture',
        plotdir='reports/figures/benchmark/'
    output:
        json='data/interim/fits/benchmark/isolated_fits.json'
    script:
        'src/fit_louvain.py'

rule benchmark_icat:
    input:
        X='data/processed/benchmark/X.csv',
        obs='data/processed/benchmark/obs.csv',
        var='data/processed/benchmark/var.csv',
        json='data/interim/fits/benchmark/isolated_fits.json',
        ncfs='data/external/benchmark_ncfs_params.json'
    params:
        treatment='benchmark',
        control='sc_celseq2',
        label='mixture',
        outdir='data/results/benchmark/icat/'
    output:
        X=protected('data/results/benchmark/icat/X.csv'),
        obs=protected('data/results/benchmark/icat/obs.csv'),
        var=protected('data/results/benchmark/icat/var.csv')
    script:
        'src/benchmark_icat.py'

rule benchmark_seurat233:
    input:
        X='data/raw/benchmark/combined/X.csv',
        obs='data/raw/benchmark/combined/obs.csv',
        json='data/interim/fits/benchmark/isolated_fits.json'
    params:
        treatment='benchmark',
        control='sc_celseq2',
        label='mixture',
        seurat=config['libraries']['seurat2.3.3']
    output:
        X='data/results/benchmark/seurat233/X.csv',
        obs='data/results/benchmark/seurat233/obs.csv'
    script:
        'src/run_seurat2-3-3.R'

rule benchmark_seurat311:
    input:
        X='data/raw/benchmark/combined/X.csv',
        obs='data/raw/benchmark/combined/obs.csv',
        json='data/interim/fits/benchmark/isolated_fits.json'
    params:
        treatment='benchmark',
        seurat=config['libraries']['seurat3.1.1']
    output:
        X='data/results/benchmark/seurat311/X.csv',
        obs='data/results/benchmark/seurat311/obs.csv',
    script:
        "src/run_seurat3-1-1.R"

rule benchmark_seurat_icat:
    input:
        X='data/results/benchmark/seurat311/X.csv',
        obs='data/results/benchmark/seurat311/obs.csv',
        json='data/interim/fits/benchmark/isolated_fits.json',
        ncfs='data/external/benchmark_ncfs_params.json'
    params:
        treatment='benchmark',
        controls='sc_celseq2',
        outdir='data/results/benchmark/seurat_icat/',
        cluster='seurat.sslouvain'
    output:
        X='data/results/benchmark/seurat_icat/X.csv',
        obs='data/results/benchmark/seurat_icat/obs.csv'
    script:
        'src/scanorama_icat.py'

rule benchmark_scanorama:
    input:
        X='data/processed/benchmark/X.csv',
        obs='data/processed/benchmark/obs.csv',
        var='data/processed/benchmark/var.csv',
        json='data/interim/fits/benchmark/isolated_fits.json'
    params:
        treatment='benchmark',
        controls='sc_celseq2',
        outdir='data/results/benchmark/scanorama/'
    output:
        X=protected('data/results/benchmark/scanorama/X.csv'),
        obs=protected('data/results/benchmark/scanorama/obs.csv'),
        var=protected('data/results/benchmark/scanorama/var.csv')
    script:
        'src/run_scanorama.py'
        
rule benchmark_scanorama_icat:
    input:
        X='data/results/benchmark/scanorama/X.csv',
        obs='data/results/benchmark/scanorama/obs.csv',
        var='data/results/benchmark/scanorama/var.csv',
        json='data/interim/fits/benchmark/isolated_fits.json',
        ncfs='data/external/benchmark_ncfs_params.json'
    output:
        X=protected('data/results/benchmark/icat_scan/X.csv'),
        obs=protected('data/results/benchmark/icat_scan/obs.csv'),
        var=protected('data/results/benchmark/icat_scan/var.csv'),
    params:
        outdir='data/results/benchmark/icat_scan/',
        treatment='benchmark',
        controls='sc_celseq2',
        cluster='scanorama.sslouvain'
    script:
        'src/scanorama_icat.py'

rule summarize_benchmark:
    input:
        obs=['data/results/benchmark/{method}/obs.csv'.format(
             method=method) for method in METHODS],
    params:
        identity='mixture'
    output:
        csv='data/results/benchmark/results.csv',
        svg='reports/figures/benchmark/metrics.svg'
    script:
        'src/summarize_benchmark.py'

# create plots/figures
rule plot_benchmark:
    input:
        results='data/results/benchmark/results.csv',
        labels=['data/results/benchmark/icat/obs.csv',
                'data/results/benchmark/seurat233/clustered.csv',
                'data/results/benchmark/scanorama/obs.csv',
                'data/results/benchmark/icat_scan/obs.csv'],
        Xs=['data/results/benchmark/icat/X.csv',
            'data/results/benchmark/icat/X.csv', # not going to use the X data
            'data/results/benchmark/scanorama/X.csv',
            'data/results/benchmark/icat_scan/X.csv'],
        fit='data/interim/fits/benchmark/isolated_fits.json'
    params:
        methods=['icat', 'seurat233', 'scanorama', 'icat_scan'],
        plotdir='reports/figures/benchmark/',
        label='mixture',
        treatment='benchmark',
        xlabel='Cell Mixture'
    output:
        metrics='reports/figures/benchmark/metrics.svg',
        method_plots=['reports/figures/benchmark/{method}/{plot}'.format(
                      method=method, plot=plot) for method, plot in\
                      itertools.product(['icat', 'seurat233', 'scanorama',
                                         'icat_scan'],
                                        ['known_cells_umap.svg',
                                         'cluster_umap.svg',
                                         'treatment_umap.svg',
                                         'cell_type_distribution.svg',
                                         'cluster_distribution.svg'])]
    script:
        'src/plot_performance.py'

# ---------------------------- Analyze Kang Data -------------------------------

# normalize counts, find most variable genes
rule kang_filter:
    input:
        X='data/processed/Kang/X.csv',
        obs='data/processed/Kang/obs.csv',
        var='data/processed/Kang/var.csv'
    output:
        X='data/processed/Kang/filtered/X.csv',
        obs='data/processed/Kang/filtered/obs.csv',
        var='data/processed/Kang/filtered/var.csv',
    params:
        outdir='data/filtered/Kang',
        plotdir='figures/Kang/'
    script:
        "src/filter_kang.py"

    
