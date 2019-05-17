

rule all:
    input:
        'data/processed/kang_2018_umi.pkl'

rule make_count_matrix:
    params:
        'data/raw'
    output:
        'data/processed/kang_2018_umi.pkl'
    script:
        'src/generate_kang_et_al.py'