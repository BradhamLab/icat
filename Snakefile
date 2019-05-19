

rule all:
    input:
        'data/processed/simulations.csv'

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
        csv='data/processed/simulations.csv'
    params:
        sims=3,
        reps=3,
        outdir='data/processed/'
    script:
        'src/generate_simulated_datasets.py'