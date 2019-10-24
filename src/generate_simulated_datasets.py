import json
import os

import numpy as np
import pandas as pd

from icat.src import simulate
from icat.src import utils
from downstream.src.analysis import utils as dutils

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
        print(sims, reps)
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
        perturbation_keys = []
        for k, v in perturbations.items():
            simmed, controls = run_simulation(control_params, v, sims, reps,
                                              controls=controls)
            exp_key = "{}.Perturbation{}".format(exp, k)
            perturbation_keys.append(exp_key)
            # keep experiment keys to combine perturbations if need be
            datasets[exp_key] = simmed
            flattened =  dict(utils.flatten_dict(control_params),
                              **utils.flatten_dict({'perturbation': v}))
            for k, v in flattened.items():
                if isinstance(v, list):
                    v = ';'.join([str(x) for x in v])
                flattened[k] = v
            flattened['dropout'] = np.sum(controls.X == 0) / controls.X.size
            csv_dict[exp_key] = flattened  # TODO: look at this with change
        # combine dataset across perturbations 
        if params['combine_perturbations']:
            # change experiment key from Experiment{Y}.Perturbation{X} to
            # Experiment{Y} given joining of perturbations
            exp_key = perturbation_keys[0].split('.')[0]
            sim_rep_data = [[[None] for i in range(reps)] for j in range(sims)]
            for sim in range(sims):
                for rep in range(reps):
                    # combine perturbations across the same simulations + reps
                    adatas = []
                    # include controls
                    pert = perturbation_keys[0]
                    adata = datasets[pert][sim][rep]
                    # rename treatment from Perturbed to Perturbed{X}
                    label = pert.split('.')[-1]
                    adata.obs['Treatment'].replace('Perturbed', label,
                                                    inplace=True)
                    adatas.append(adata)
                    for pert in perturbation_keys[1:]:
                        adata = datasets[pert][sim][rep]
                        adata = dutils.filter_cells(adata, 'Treatment',
                                                    lambda x: x == 'Perturbed')
                        # rename treatment from Perturbed to Perturbed{X}
                        label = pert.split('.')[-1]
                        adata.obs['Treatment'].replace('Perturbed', label,
                                                       inplace=True)
                        adatas.append(adata)
                    adata = utils.rbind_adata(adatas)
                    adata.obs.index = ['cell-{}'.format(i + 1)\
                                       for i in range(adata.shape[0])]
                    sim_rep_data[sim][rep] = adata
            # remove uncombined data from dataset dict, add combined with new key
            for each in perturbation_keys:
                datasets.pop(each)
                datasets[exp_key] = sim_rep_data

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
        for n_sim in range(sims):
            replicates = datasets[each][n_sim]
            for n_rep in range(reps):
                data = replicates[n_rep]
                exp_dir = os.path.join(out_dir,
                                       "{}Sim{}Rep{}".format(each, n_sim + 1,
                                                             n_rep + 1))
                if not os.path.exists(exp_dir):
                    os.makedirs(exp_dir)
                # do not write over previously created data
                if not os.path.exists(os.path.join(exp_dir, 'X.csv')):
                    data.write_csvs(dirname=exp_dir, skip_data=False)
                    
