
import glob
import os
import itertools
from icat.src import snakemake_utils as utils

configfile: './snakemake-config.yaml'

FILES = ['X.csv', 'obs.csv', 'var.csv']
CONDITIONS = ['Controls', 'Treated']
EXPERIMENTS = utils.get_simulation_ids(config['simulations']['json'],
                                       config['simulations']['sims'],
                                       config['simulations']['reps'])

SIMULATED = ["data/processed/simulated/{exp}/{treat}/{out}".\
             format(exp=exp, treat=treat, out=out)\
             for exp, treat, out in itertools.product(EXPERIMENTS,
                                                      CONDITIONS, FILES)]
BENCHMARK = ['cellmix1', 'cellmix2', 'cellmix3', 'cellmix4', 'sc_celseq2']
MIXES = BENCHMARK[:-1]

rule all:
    input:
        ['data/results/clustered/icat/{exp}/performance.csv'.format(exp=exp)\
            for exp in EXPERIMENTS],
        ['data/results/clustered/seurat/{exp}/clustered.csv'.format(exp=exp)\
            for exp in EXPERIMENTS],
        ['data/results/clustered/scanorama/{exp}/obs.csv'.format(exp=exp)\
            for exp in EXPERIMENTS],
        'data/processed/Kang/X.csv'
        # ['data/processed/BenchData/{bench}/X.csv'.format(bench=bench)\
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

rule format_benchmark_data:
    input:
        counts=['data/raw/BenchData/{bench}.count.csv'.format(bench=bench)\
                for bench in BENCHMARK],
        meta=['data/raw/BenchData/{bench}.metadata.csv'.format(bench=bench)\
              for bench in BENCHMARK]
    params:
        outdir='data/processed/BenchData/',
        plotdir='figures/benchmark/'
    output:
        X=['data/processed/BenchData/{bench}/X.csv'.format(bench=bench)
           for bench in BENCHMARK],
        obs=['data/processed/BenchData/{bench}/obs.csv'.format(bench=bench)
             for bench in BENCHMARK],
        gene_svg='figures/benchmark/filter_genes_dispersion.svg'
    script:
        'src/format_benchmark_data.py'

rule concatenate_benchmark_mixtures:
    input:
        X=['data/processed/BenchData/{mix}/X.csv'.format(mix=mix)
           for mix in MIXES],
        obs=['data/processed/BenchData/{mix}/obs.csv'.format(mix=mix)
             for mix in MIXES],
        var=['data/processed/BenchData/{mix}/var.csv'.format(mix=mix)
             for mix in MIXES]
    params:
        outdir='data/processed/BenchData/mixed/'
    output:
        X='data/processed/BenchData/mixed/X.csv',
        obs='data/processed/BenchData/mixed/obs.csv',
        var='data/processed/BenchData/mixed/var.csv'
    script:
        "src/concatenate_mixtures.py"

rule move_benchmark_isolated:
    input:
        X='data/processed/BenchData/sc_celseq2/X.csv',
        obs='data/processed/BenchData/sc_celseq2/obs.csv',
        var='data/processed/BenchData/sc_celseq2/var.csv'
    output:
        X='data/processed/BenchData/isolated/X.csv',
        obs='data/processed/BenchData/isolated/obs.csv',
        var='data/processed/BenchData/isolated/var.csv'
    shell:
        """
        mv {input.X} {output.X};
        mv {input.obs} {output.obs};
        mv {input.var} {output.var} 
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
        ctrl_X='data/processed/simulated/{exp}/Controls/X.csv',
        ctrl_obs='data/processed/simulated/{exp}/Controls/obs.csv',
        ctrl_var='data/processed/simulated/{exp}/Controls/var.csv',
        prtb_X='data/processed/simulated/{exp}/Treated/X.csv',
        prtb_obs='data/processed/simulated/{exp}/Treated/obs.csv',
        prtb_var='data/processed/simulated/{exp}/Treated/var.csv'
    params:
        label='Population',
        plotdir='figures/simulated/{exp}/'
    output:
        json='data/interim/fits/{exp}Controls_fits.json',
        ctrl_svg='figures/simulated/{exp}/umap_controls.svg',
        prtb_svg='figures/simulated/{exp}/umap_treated.svg',
        comb_svg='figures/simulated/{exp}/umap_combined.svg'
    script:
        'src/fit_louvain.py'

rule simulated_icat:
    input:
        ctrl_X='data/processed/simulated/{exp}/Controls/X.csv',
        ctrl_obs='data/processed/simulated/{exp}/Controls/obs.csv',
        ctrl_var='data/processed/simulated/{exp}/Controls/var.csv',
        prtb_X='data/processed/simulated/{exp}/Treated/X.csv',
        prtb_obs='data/processed/simulated/{exp}/Treated/obs.csv',
        prtb_var='data/processed/simulated/{exp}/Treated/var.csv',
        json='data/interim/fits/{exp}Controls_fits.json'
    params:
        name='{exp}',
        plotdir='figures/clustered/{exp}/',
        outdir='data/results/clustered/icat/{exp}',
    output:
        csv='data/results/clustered/icat/{exp}/performance.csv',
        obs='data/results/clustered/icat/{exp}/obs.csv',
        p1='figures/clustered/{exp}/umap_louvain.svg',
        p2='figures/clustered/{exp}/umap_ncfs-louvain.svg',
        p3='figures/clustered/{exp}/umap_sslouvain.svg'
    script:
        'src/evaluate_icat.py'

rule simulated_seurat:
    input:
        ctrl_X='data/processed/simulated/{exp}/Controls/X.csv',
        ctrl_obs='data/processed/simulated/{exp}/Controls/obs.csv',
        prtb_X='data/processed/simulated/{exp}/Treated/X.csv',
        prtb_obs='data/processed/simulated/{exp}/Treated/obs.csv',
        json='data/interim/fits/{exp}Controls_fits.json'
    params:
        name='{exp}'
    output:
        csv='data/results/clustered/seurat/{exp}/clustered.csv'
    script:
        'src/cluster_seurat.R'

rule simulated_scanorama:
    input:
        ctrl_X='data/processed/simulated/{exp}/Controls/X.csv',
        ctrl_obs='data/processed/simulated/{exp}/Controls/obs.csv',
        ctrl_var='data/processed/simulated/{exp}/Controls/var.csv',
        prtb_X='data/processed/simulated/{exp}/Treated/X.csv',
        prtb_obs='data/processed/simulated/{exp}/Treated/obs.csv',
        prtb_var='data/processed/simulated/{exp}/Treated/var.csv',
        json='data/interim/fits/{exp}Controls_fits.json'
    output:
        X='data/results/clustered/scanorama/{exp}/X.csv',
        obs='data/results/clustered/scanorama/{exp}/obs.csv',
        var='data/results/clustered/scanorama/{exp}/var.csv'
    params:
        name='{exp}',
        outdir='data/results/clustered/scanorama/{exp}/'
    script:
        'src/evaluate_scanorama.py'
    
# ------------------------ Analyze Benchmark Data ------------------------------

rule fit_benchmark_data:
    input:
        ctrl_X='data/processed/BenchData/sc_celseq2/X.csv',
        ctrl_obs='data/processed/BenchData/sc_celseq2/obs.csv',
        prtb_X = ['data/processed/BenchData/{}/X.csv'.format(x)\
                  for x in ['cellmix1', 'cellmix2', 'cellmix3', 'cellmix4']],
        prtb_obs = ['data/processed/BenchData/{}/obs.csv'.format(x)\
                    for x in ['cellmix1', 'cellmix2', 'cellmix3', 'cellmix4']]
    params:
        label='cell_line',
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
        X=['data/processed/BenchData/{bench}/X.csv'.format(bench=bench)\
           for bench in BENCHMARK],
        obs=['data/processed/BenchData/{bench}/obs.csv'.format(bench=bench)\
             for bench in BENCHMARK],
        json='data/interim/fits/benchmark/isolated_fits.json'
    params:
        control_id='sc_celseq2'
    output:
        'data/results/benchmark/icat_clusters.csv'
    script:
        'src/benchmark_icat.py'

rule benchmark_seurat:
    input:
        ctrl_X='data/processed/BenchData/isolated/X.csv',
        ctrl_obs='data/processed/BenchData/isolated/obs.csv',
        prtb_X='data/processed/BenchData/mixed/X.csv',
        prtb_obs='data/processed/BenchData/mixed/obs.csv',
        json='data/interim/fits/benchmark/isolated_fits.json'
    params:
        name='benchmark'
    output:
        csv='data/results/clustered/seurat/benchmark/clustered.csv'
    script:
        'src/cluster_seurat.R'

rule benchmark_scanorama:
    input:
        ctrl_X='data/processed/BenchData/isolated/X.csv',
        ctrl_obs='data/processed/BenchData/isolated/obs.csv',
        ctrl_var='data/processed/Benchdata/isolated/Controls/var.csv',
        ctrl_X='data/processed/BenchData/mixed/X.csv',
        ctrl_obs='data/processed/BenchData/mixed/obs.csv',
        ctrl_var='data/processed/Benchdata/mixed/Controls/var.csv',
        json='data/interim/fits/benchmark/isolated_fits.json'
    output:
        X='data/results/clustered/scanorama/benchmark/X.csv',
        obs='data/results/clustered/scanorama/benchmark/obs.csv',
        var='data/results/clustered/scanorama/benchmark/var.csv'
    params:
        name='benchmark',
        outdir='data/results/clustered/scanorama/benchmark/'
    script:
        'src/evaluate_scanorama.py'


rule benchmark_scanorama_icat:
    input:
        X='data/results/clustered/scanorama/benchmark/X.csv',
        obs='data/results/clustered/scanorama/benchmark/obs.csv',
        var='data/results/clustered/scanorama/benchmark/var.csv',
        json=json='data/interim/fits/benchmark/isolated_fits.json'
    output:
        'data/results/clustered/icat_scan/benchmark/'
    params:
        outdir:'data/results/clustered/icat_scan/benchmark/',
        treat_col:'benchmark',
        treat_values = MIXES,
    script:
        'src/scanorama_icat.py'
    
