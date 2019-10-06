
import glob
import os
import itertools
from icat.src import snakemake_utils as utils

configfile: './snakemake-config.yaml'

FILES = ['X.csv', 'obs.csv', 'var.csv']
EXPERIMENTS = utils.get_simulation_ids(config['simulations']['json'],
                                       config['simulations']['sims'],
                                       config['simulations']['reps'])

SIMULATED = ["data/processed/simulated/{exp}/{out}".\
             format(exp=exp, out=out)\
             for exp, out in itertools.product(EXPERIMENTS, FILES)]
BENCHMARK = ['cellmix1', 'cellmix2', 'cellmix3', 'cellmix4', 'sc_celseq2']
MIXES = BENCHMARK[:-1]

rule all:
    input:
        ['data/results/simulated/icat/{exp}/performance.csv'.format(exp=exp)\
            for exp in EXPERIMENTS],
        ['data/results/simulated/seurat/{exp}/clustered.csv'.format(exp=exp)\
            for exp in EXPERIMENTS],
        ['data/results/simulated/scanorama/{exp}/obs.csv'.format(exp=exp)\
            for exp in EXPERIMENTS],
        'data/processed/Kang/X.csv',
        'data/results/benchmark/results.csv'
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
        X=['data/processed/benchmark/{bench}/X.csv'.format(bench=bench)
           for bench in BENCHMARK],
        obs=['data/processed/benchmark/{bench}/obs.csv'.format(bench=bench)
             for bench in BENCHMARK],
        var=['data/processed/benchmark/{bench}/var.csv'.format(bench=bench)
             for bench in BENCHMARK], 
        gene_svg='figures/benchmark/filter_genes_dispersion.svg'
    script:
        'src/format_benchmark_data.py'

rule concatenate_benchmark_mixtures:
    input:
        X=['data/processed/benchmark/{mix}/X.csv'.format(mix=mix)
           for mix in MIXES],
        obs=['data/processed/benchmark/{mix}/obs.csv'.format(mix=mix)
             for mix in MIXES],
        var=['data/processed/benchmark/{mix}/var.csv'.format(mix=mix)
             for mix in MIXES]
    params:
        outdir='data/processed/benchmark/mixed/'
    output:
        X='data/processed/benchmark/mixed/X.csv',
        obs='data/processed/benchmark/mixed/obs.csv',
        var='data/processed/benchmark/mixed/var.csv'
    script:
        "src/concatenate_mixtures.py"

rule move_benchmark_isolated:
    input:
        X='data/processed/benchmark/sc_celseq2/X.csv',
        obs='data/processed/benchmark/sc_celseq2/obs.csv',
        var='data/processed/benchmark/sc_celseq2/var.csv'
    output:
        X='data/processed/benchmark/isolated/X.csv',
        obs='data/processed/benchmark/isolated/obs.csv',
        var='data/processed/benchmark/isolated/var.csv'
    shell:
        """
        cp {input.X} {output.X};
        cp {input.obs} {output.obs};
        cp {input.var} {output.var} 
        """
    
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
        X='data/processed/simulated/{exp}/X.csv',
        obs='data/processed/simulated/{exp}/obs.csv',
        var='data/processed/simulated/{exp}/var.csv'
    params:
        treatment='Treatment',
        control='Control',
        label='Population',
        plotdir='reports/figures/simulated/{exp}/'
    output:
        json='data/interim/fits/simulated/{exp}_fit.json'
    script:
        'src/fit_louvain.py'

rule simulated_icat:
    input:
        X='data/processed/simulated/{exp}/X.csv',
        obs='data/processed/simulated/{exp}/obs.csv',
        var='data/processed/simulated/{exp}/var.csv',
        json='data/interim/fits/simulated/{exp}_fit.json',
        ncfs='data/external/simulated_ncfs_params.json'
    params:
        name='{exp}',
        treatment='Treatment',
        control='Control',
        plotdir='reports/figures/simulated/{exp}/icat/',
        outdir='data/results/simulated/icat/{exp}/',
    output:
        csv='data/results/simulated/icat/{exp}/performance.csv',
        obs='data/results/simulated/icat/{exp}/obs.csv',
        p1='reports/figures/simulated/{exp}/icat/umap_Louvain.svg',
        p2='reports/figures/simulated/{exp}/icat/umap_NCFS-Louvain.svg',
        p3='reports/figures/simulated/{exp}/icat/umap_NCFS-SSLouvain.svg'
    script:
        'src/evaluate_icat.py'

rule simulated_seurat:
    input:
        X='data/processed/simulated/{exp}/X.csv',
        obs='data/processed/simulated/{exp}/obs.csv',
        json='data/interim/fits/simulated/{exp}_fit.json'
    params:
        treatment='Treatment',
        control='Control',
        label='Population',
        seurat=config['libraries']['seurat2.3.3']
    output:
        csv='data/results/simulated/seurat/{exp}/clustered.csv'
    script:
        'src/cluster_seurat.R'

rule simulated_scanorama:
    input:
        X='data/processed/simulated/{exp}/X.csv',
        obs='data/processed/simulated/{exp}/obs.csv',
        var='data/processed/simulated/{exp}/var.csv',
        json='data/interim/fits/simulated/{exp}_fit.json'
    output:
        X='data/results/simulated/scanorama/{exp}/X.csv',
        obs='data/results/simulated/scanorama/{exp}/obs.csv',
        var='data/results/simulated/scanorama/{exp}/var.csv'
    params:
        treatment='Treatment',
        controls='Control',
        outdir='data/results/simulated/scanorama/{exp}/'
    script:
        'src/run_scanorama.py'

rule simulated_scanorama_icat:
    input:
        X='data/results/simulated/scanorama/{exp}/X.csv',
        obs='data/results/simulated/scanorama/{exp}/obs.csv',
        var='data/results/simulated/scanorama/{exp}/var.csv',
        json='data/interim/fits/simulated/{exp}Controls_fit.json',
        ncfs='data/external/benchmark_ncfs_params.json'
    output:
        X='data/results/simulated/icat_scan/{exp}/X.csv',
        obs='data/results/simulated/icat_scan/{exp}/obs.csv',
        var='data/results/simulated/icat_scan/{exp}/var.csv'
    params:
        outdir='data/results/benchmark/icat_scan/{exp}/',
        treat_col='Treatment',
        treat_values = MIXES
    script:
        'src/scanorama_icat.py'

rule summarize_simulated:
    input:
        perf='data/results/simulated/icat/{exp}/performance.csv',
        icat='data/results/simulated/icat/{exp}/obs.csv',
        gene='data/results/simulated/icat/{exp}/var.csv',
        seurat='data/results/simulated/seurat/{exp}/clustered.csv',
        scanorama='data/results/simulated/scanorama/{exp}/obs.csv',
    output:
        performance='data/results/simulated/performance.csv'
    params:
        plotdir='reports/figures/simulated/{exp}/'
    script:
        'src/summarize_simulated.py'
    
# ------------------------ Analyze Benchmark Data ------------------------------
rule fit_benchmark_data:
    input:
        ctrl_X='data/processed/benchmark/sc_celseq2/X.csv',
        ctrl_obs='data/processed/benchmark/sc_celseq2/obs.csv',
        prtb_X = ['data/processed/benchmark/{}/X.csv'.format(x)\
                  for x in ['cellmix1', 'cellmix2', 'cellmix3', 'cellmix4']],
        prtb_obs = ['data/processed/benchmark/{}/obs.csv'.format(x)\
                    for x in ['cellmix1', 'cellmix2', 'cellmix3', 'cellmix4']]
    params:
        label='mixture',
        plotdir='figures/benchmark/'
    output:
        json='data/interim/fits/benchmark/isolated_fits.json',
        ctrl_svg='figures/benchmark/umap_isolated_cells.svg',
        prtb_svg='figures/benchmark/umap_mixed_cells.svg',
        comb_svg='figures/benchmark/umap_combined.svg'
    script:
        'src/fit_louvain.py'

rule benchmark_icat:
    input:
        X=['data/processed/benchmark/{bench}/X.csv'.format(bench=bench)\
           for bench in BENCHMARK],
        obs=['data/processed/benchmark/{bench}/obs.csv'.format(bench=bench)\
             for bench in BENCHMARK],
        json='data/interim/fits/benchmark/isolated_fits.json',
        ncfs='data/external/benchmark_ncfs_params.json'
    params:
        control_id='sc_celseq2',
        outdir='data/results/benchmark/icat/'
    output:
        X='data/results/benchmark/icat/X.csv',
        obs='data/results/benchmark/icat/obs.csv',
        var='data/results/benchmark/icat/var.csv'
    script:
        'src/benchmark_icat.py'

rule benchmark_seurat:
    input:
        ctrl_X='data/processed/benchmark/isolated/X.csv',
        ctrl_obs='data/processed/benchmark/isolated/obs.csv',
        prtb_X='data/processed/benchmark/mixed/X.csv',
        prtb_obs='data/processed/benchmark/mixed/obs.csv',
        json='data/interim/fits/benchmark/isolated_fits.json'
    params:
        name='benchmark'
    output:
        csv='data/results/benchmark/seurat/clustered.csv'
    script:
        'src/cluster_seurat.R'

rule benchmark_scanorama:
    input:
        X=['data/processed/benchmark/{bench}/X.csv'.format(bench=bench)\
           for bench in BENCHMARK],
        obs=['data/processed/benchmark/{bench}/obs.csv'.format(bench=bench)\
             for bench in BENCHMARK],
        var=['data/processed/benchmark/{bench}/var.csv'.format(bench=bench)\
             for bench in BENCHMARK],
        json='data/interim/fits/benchmark/isolated_fits.json'
    params:
        control_id='sc_celseq2',
        outdir='data/results/benchmark/scanorama/'
    output:
        X='data/results/benchmark/scanorama/X.csv',
        obs='data/results/benchmark/scanorama/obs.csv',
        var='data/results/benchmark/scanorama/var.csv'
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
        X='data/results/benchmark/icat_scan/X.csv',
        obs='data/results/benchmark/icat_scan/obs.csv',
        var='data/results/benchmark/icat_scan/var.csv'
    params:
        outdir='data/results/benchmark/icat_scan/',
        treat_col='benchmark',
        treat_values = MIXES
    script:
        'src/scanorama_icat.py'

rule summarize_benchmark:
    input:
        icat='data/results/benchmark/icat/obs.csv',
        seurat='data/results/benchmark/seurat/clustered.csv',
        scanorama='data/results/benchmark/scanorama/obs.csv',
        icat_scan='data/results/benchmark/icat_scan/obs.csv',
    params:
        identity='mixture'
    output:
        csv='data/results/benchmark/results.csv'
    script:
        'src/summarize_benchmark.py'


# create plots/figures
rule plot_benchmark:
    input:
        results='data/results/benchmark/results.csv',
        labels=['data/results/benchmark/icat/obs.csv',
                'data/results/benchmark/seurat/clustered.csv',
                'data/results/benchmark/scanorama/obs.csv',
                'data/results/benchmark/icat_scan/obs.csv'],
        Xs=['data/results/benchmark/icat/X.csv',
            'data/results/benchmark/icat/X.csv', # not going to use the X data
            'data/results/benchmark/scanorama/X.csv',
            'data/results/benchmark/icat_scan/X.csv'],
        fit='data/interim/fits/benchmark/isolated_fits.json'
    params:
        methods=['icat', 'seurat', 'scanorama', 'icat_scan'],
        plotdir='reports/figures/benchmark/',
        label='mixture',
        treatment='benchmark',
        xlabel='Cell Mixture'
    output:
        metrics='reports/figures/benchmark/metrics.svg',
        method_plots=['reports/figures/benchmark/{method}/{plot}'.format(
                      method=method, plot=plot) for method, plot in\
                      itertools.product(['icat', 'seurat', 'scanorama',
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
        ctrl_X='data/filtered/Kang/controls/X.csv',
        ctrl_obs='data/filtered/Kang/controls/obs.csv',
        ctrl_var='data/filtered/Kang/controls/var.csv',
        treated_X='data/filtered/Kang/treated/X.csv',
        treated_obs='data/filtered/Kang/treated/obs.csv',
        treated_var='data/filtered/Kang/treated/var.csv',
    params:
        outdir='data/filtered/Kang',
        plotdir='figures/Kang/'
    script:
        "src/filter_kang.py"

    
