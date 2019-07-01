
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
        ['data/results/{exp}_icat_performance.csv'.format(exp=exp)\
        for exp in EXPERIMENTS]

rule simulate_data:
    input:
        input_json=os.path.join('../../', config['simulations']['json']),
    params:
        sims=config['simulations']['sims'],
        reps=config['simulations']['reps'],
        outdir="data/processed/simulated/"
    output:
        data=SIMULATED,
        csv="../../data/processed/simulated/simulations.csv"
    script:
        "../../src/generate_simulated_datasets.py"

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
        plotdir='figures/clustered/{exp}/'
    output:
        csv='data/results/{exp}_icat_performance.csv',
        p1='figures/clustered/{exp}/umap_louvain.svg',
        p2='figures/clustered/{exp}/umap_ncfs-louvain.svg',
        p3='figures/clustered/{exp}/umap_sslouvain.svg'
    script:
        'src/evaluation.py'


rule cluster_seurat:
    input:
        ctrl_X='data/processed/simulated/{exp}/Controls/X.csv',
        ctrl_obs='data/processed/simulated/{exp}/Controls/obs.csv',
        ctrl_var='data/processed/simulated/{exp}/Controls/var.csv',
        prtb_X='data/processed/simulated/{exp}/Treated/X.csv',
        prtb_obs='data/processed/simulated/{exp}/Treated/obs.csv',
        prtb_var='data/processed/simulated/{exp}/Treated/var.csv',
        json='data/interim/fits/{exp}Controls_fits.json'
    output:
        csv='data/interim/seurat/{exp}/clustered.csv'
    params:
        name='{sim}'
    script:
        'src/cluster_seurat.R'
    