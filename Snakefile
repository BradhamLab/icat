

rule all:
    input:
        'data/interim/control_louvain_fits.csv'

rule make_count_matrix:
    params:
        'data/raw'
    output:
        'data/processed/kang_2018_umi.pkl'
    script:
        'src/generate_kang_et_al.py'

rule simulate_data:
    input:
        json='data/external/experiments.json'
    output:
        csv='data/processed/simulated/simulations.csv',
    params:
        sims=3,
        reps=3,
        outdir='data/processed/simulated/'
    script:
        'src/generate_simulated_datasets.py'

rule fit_louvain:
    input:
        'data/processed/simulated/simulations.csv'
    params:
        data = 'data/processed/simulated',
        label = 'Population',
        regex = '*Controls.pkl',
        plotdir = 'plots/simulated/'
    output:
        csv = 'data/interim/control_louvain_fits.csv'
    script:
        'src/fit_louvain.py'