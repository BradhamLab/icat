import json
import os
import pickle as pkl

import pandas as pd

from icat.src import simulate
from icat.src import utils

def parse_params(sim_params):
    if 'dispersion' in sim_params and sim_params['dispersion'] == 'random':
        sim_params['dispersion'] = simulate.dispersions(sim_params['genes'])
    if 'gene_targets' in sim_params:
        if sim_params['gene_targets'] == 'None':
            sim_params['gene_targets'] = None
    return sim_params

def run_simulation(c_params, p_params, sims, reps, controls=None):
        c_params = parse_params(c_params)
        p_params = parse_params(p_params)
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
        return (experiment.run(simulations=sims, replications=reps,
                               controls=controls, pop_targets=pop_targets),
                controls)

def main(configs, sims=1, reps=1):
    datasets = dict()
    csv_dict = dict()
    for exp, params in configs.items():
        control_params = params['controls']
        perturbations = params['perturbations']
        controls = None
        for k, v in perturbations.items():
            simmed, controls = run_simulation(control_params, v, sims, reps,
                                              controls=controls)
            exp_key = "{}.Perturbation{}".format(exp, k)
            datasets[exp_key] = simmed
            flattened =  dict(utils.flatten_dict(control_params),
                              **utils.flatten_dict({'perturbation': v}))
            for k, v in flattened.items():
                if isinstance(v, list):
                    v = ';'.join([str(x) for x in v])
                flattened[k] = v
            csv_dict[exp_key] = flattened  # TODO: look at this with change
        # parameters defining control space

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
                exp_dir = os.path.join(out_dir,
                                       "{}Sim{}Rep{}".format(each, n_sim + 1,
                                                             n_rep + 1))
                for i, x in enumerate(['Controls', 'Treated']):
                    write_dir = os.path.join(exp_dir, x)
                    data[i].write_csvs(dirname=write_dir, skip_data=False)
                    
