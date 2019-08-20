
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

rule all:
    input:
        ['data/results/clustered/icat/{exp}/performance.csv'.format(exp=exp)\
            for exp in EXPERIMENTS],
        ['data/results/clustered/seurat/{exp}/clustered.csv'.format(exp=exp)\
            for exp in EXPERIMENTS],
        ['data/results/clustered/scanorama/{exp}/obs.csv'.format(exp=exp)\
            for exp in EXPERIMENTS],

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

rule fit_louvain:
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

rule cluster_icat:
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

rule cluster_seurat:
    input:
        ctrl_X='data/processed/simulated/{exp}/Controls/X.csv',
        ctrl_obs='data/processed/simulated/{exp}/Controls/obs.csv',
        ctrl_var='data/processed/simulated/{exp}/Controls/var.csv',
        prtb_X='data/processed/simulated/{exp}/Treated/X.csv',
        prtb_obs='data/processed/simulated/{exp}/Treated/obs.csv',
        prtb_var='data/processed/simulated/{exp}/Treated/var.csv',
        json='data/interim/fits/{exp}Controls_fits.json'
    params:
        name='{exp}'
    output:
        csv='data/results/clustered/seurat/{exp}/clustered.csv'
    script:
        'src/cluster_seurat.R'

rule cluster_scanorama:
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
    