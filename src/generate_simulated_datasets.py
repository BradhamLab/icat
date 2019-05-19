import json
import os
import pickle as pkl

import pandas as pd

import simulate
import utils


def main(configs, sims=1, reps=1):
    datasets = dict()
    csv_dict = dict()
    for exp, params in configs.items():
        values = utils.flatten_dict(params)
        for k, v in values.items():
            if isinstance(v, list):
                v = ';'.join([str(x) for x in v])
            values[k] = v
        csv_dict[exp] = values
        c_params = params['controls']
        p_params = params['perturbation']
        if c_params['dispersion'] == 'random':
            c_params['dispersion'] = simulate.dispersions(c_params['genes'])
        # not interested in the current functionality of pop_targets
        pop_targets = p_params.pop('pop_targets')
        if pop_targets == 'None':
            pop_targets = None
        if p_params['gene_targets'] == 'None':
            p_params['gene_targets'] = None
        experiment = simulate.Experiment(control_kwargs=c_params,
                                         perturb_kwargs=p_params)
        datasets[exp] = experiment.run(simulations=sims, replications=reps,
                                       pop_targets=pop_targets)
    return datasets, pd.DataFrame(csv_dict).T
    

if __name__ == '__main__':
    input_json = '../data/external/experiments.json'
    sims = 1
    reps = 1
    out_csv = '../data/raw/simulation_paramters.csv'
    out_dir = '../data/raw/'
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        input_json = snakemake.input['json']
        sims = snakemake.params['sims']
        reps = snakemake.params['reps']
        out_csv = snakemake.output['csv']
        out_dir = snakemake.params['outdir']
    with open(input_json) as f:
        configs = json.load(f)
    datasets, sim_data = main(configs, sims, reps)
    sim_data.to_csv(out_csv)
    for each in datasets.keys():
        for n_sim in range(len(datasets[each])):
            independents = datasets[each][n_sim]
            for n_rep in range(len(independents)):
                data = [independents[n_rep]['controls'],
                        independents[n_rep]['treated']]
                for i, x in enumerate(['Controls.pkl', 'Treated.pkl']):
                    fn = "{}Sim{}Rep{}-{}".format(each, n_sim + 1, n_rep + 1, x)
                    with open(os.path.join(out_dir, fn), 'wb') as f:
                        pkl.dump(data[i], f)
