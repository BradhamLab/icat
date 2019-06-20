import json
import os
import pickle as pkl

import pandas as pd

import simulate
import utils

def parse_params(sim_params):
    if sim_params['dispersion'] == 'random':
        sim_params['dispersion'] = simulate.dispersions(sim_params['genes'])
    if 'gene_targets' in sim_params:
        if sim_params['gene_targets'] == 'None':
            sim_params['gene_targets'] = None
    return sim_params

def main(configs, sims=1, reps=1):
    datasets = dict()
    csv_dict = dict()
    for exp, params in configs.items():
        values = utils.flatten_dict(params)
        for k, v in values.items():
            if isinstance(v, list):
                v = ';'.join([str(x) for x in v])
            values[k] = v
        csv_dict[exp] = values  # TODO: look at this with change
        # parameters defining control space
        c_params = parse_params(params['controls'])
        # list of perturbations to apply to control cells 
        perturb_list = params['perturbations']
        controls = None
        for i, p_params in enumerate(perturb_list):
            # not interested in the current functionality of pop_targets in
            # SingleCellDataset.simulate(), pass to Experiment.run() instead.
            pop_targets = None
            if 'pop_targets' in p_params:
                pop_targets = p_params.pop('pop_targets')
                if pop_targets == 'None':
                    pop_targets = None
            # create experiment object for control-perturbation pairing
            experiment = simulate.Experiment(control_kwargs=c_params,
                                             perturb_kwargs=p_params)
            # simulate baseline control dataset to be used for all associated
            # perturbations
            if controls is None:
                controls = experiment.simulate_controls()
            exp_key = "{}.P{}".format(exp, i + 1)
            datasets[exp_key] = experiment.run(simulations=sims,
                                               replications=reps,
                                               controls=controls,
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
            replicates = datasets[each][n_sim]
            for n_rep in range(len(replicates)):
                data = [replicates[n_rep]['controls'],
                        replicates[n_rep]['treated']]
                for i, x in enumerate(['Controls.pkl', 'Treated.pkl']):
                    fn = "{}Sim{}Rep{}-{}".format(each, n_sim + 1, n_rep + 1, x)
                    with open(os.path.join(out_dir, fn), 'wb') as f:
                        pkl.dump(data[i], f)
