
import glob
import os

def id_from_data(data_dir, str_pattern):
    to_glob = os.path.join(data_dir, '*' + str_pattern)
    ids = [os.path.basename(x).replace(str_pattern, '')\
           for x in glob.glob(to_glob)]
    return ids
SIMS = id_from_data('data/processed/simulated', 'Controls.pkl')

rule all:
    input:
        ['data/results/{sim}_icat_performance.csv''.format(sim=sim)\
        for sim in SIMS]

rule fit_louvain:
    input:
        ctrl='data/processed/simulated/{sim}Controls.pkl',
        prtb='data/processed/simulated/{sim}Treated.pkl'
    params:
        label='Population',
        plotdir='figures/simulated/{sim}/'
    output:
        json='data/interim/fits/{sim}Controls_fits.json',
        ctrl_png='figures/simulated/{sim}/umap_controls.png',
        prtb_png='figures/simulated/{sim}/umap_treated.png',
        comb_png='figures/simulated/{sim}/umap_combined.png'
    script:
        'src/fit_louvain.py'

rule evaluate_icat:
    input:
        ctrl='data/processed/simulated/{sim}Controls.pkl',
        prtb='data/processed/simulated/{sim}Treated.pkl',
        json='data/interim/fits/{sim}Controls_fits.json'
    params:
        name='{sim}',
        plotdir='figures/clustered/{sim}/'
    output:
        csv='data/results/{sim}_icat_performance.csv',
        p1='figures/clustered/{sim}/umap_louvain.png',
        p2='figures/clustered/{sim}/umap_ncfs-louvain.png',
        p3='figures/clustered/{sim}/umap_sslouvain.png'
    script:
        'src/evaluation.py'
    