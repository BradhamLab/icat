
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
        ['figures/clustered/{sim}/seurat.svg'.format(sim=sim)\
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
        ctrl_svg='figures/simulated/{sim}/umap_controls.svg',
        prtb_svg='figures/simulated/{sim}/umap_treated.svg',
        comb_svg='figures/simulated/{sim}/umap_combined.svg'
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
        p1='figures/clustered/{sim}/umap_louvain.svg',
        p2='figures/clustered/{sim}/umap_ncfs-louvain.svg',
        p3='figures/clustered/{sim}/umap_sslouvain.svg'
    script:
        'src/evaluation.py'

rule export_data:
    input:
        ctrl='data/processed/simulated/{sim}Controls.pkl',
        prtb='data/processed/simulated/{sim}Treated.pkl',
    params:
        ctrl_dir='data/interim/csvs/{sim}/Controls/',
        prtb_dir='data/interim/csvs/{sim}/Treated/'
    output:
        'data/interim/csvs/{sim}/Controls/X.csv',
        'data/interim/csvs/{sim}/Controls/obs.csv',
        'data/interim/csvs/{sim}/Treated/X.csv',
        'data/interim/csvs/{sim}/Treated/obs.csv'
    script:
        'src/to_csv.py'

rule cluster_seurat:
    input:
        X_ctrl='data/interim/csvs/{sim}/Controls/X.csv',
        obs_ctrl='data/interim/csvs/{sim}/Controls/obs.csv',
        X_prtb='data/interim/csvs/{sim}/Treated/X.csv',
        obs_prtb='data/interim/csvs/{sim}/Treated/obs.csv'
    output:
        csv='data/interim/seurat/{sim}/clustered.csv'
    script:
        'src/cluster_seurat.R'

rule evaluate_seurat:
    input:
        csv='data/interim/seurat/{sim}/clustered.csv'
    output:
        csv='data/results/{sim}_seurat_performance.csv',
        svg='figures/clustered/{sim}/seurat.svg'
    params:
        name='{sim}'
    script:
        'src/evaluate_seurat.py'

rule concat_seurat:
    input:
        csvs=['data/results/{sim}_seurat_performance.csv'.format(sim=sim)\
              for sim in SIMS]
    output:
        csv='data/final/seurat_performance.csv'
    script:
        'src/concatenate_results.py'

    