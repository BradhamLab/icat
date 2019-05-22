
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
        ['data/interim/fits/{sim}Controls.csv'.format(sim=sim)\
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

# rule evaluate_icat:
#     input:
#         csv='data/interim/control_louvain_fits.csv'
#     params:
#         data = 'data/processed/simulated',
#         ctrl_str = "Controls",
#         prtb_str = 'Treated'
#     output:
#         csv='data/results/icat_performance.csv'
#     script:
#         'src/evaluation.py'
    