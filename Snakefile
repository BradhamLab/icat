

rule all:
    input:
        'data/processed/kang_2018_umi.pkl'

rule make_count_matrix:
    params:
        'data/raw'
    output:
        'data/processed/kang_2018_umi.pkl'
    script:
        'src/make_datasets.py'