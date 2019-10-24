import json
import re

def get_simulation_ids(exp_json, sims, reps):
    """
    Get expected ids for simulated datasets.
    
    Parameters
    ----------
    exp_json : str
        Json file containing specifications for simulation experiments.
    sims : int
        Number of times to simulate each experiment.
    reps : int
        Number of times to replicate each `sim` in simulation experiments.
    
    Returns
    -------
    list
        list of ids for each dataset expected to be simulated. 
    """
    experiments = []
    with open(exp_json, 'r') as f:
        exp_configs = json.load(f)
    for key, value in exp_configs.items():
        for perturbation in value['perturbations']:
            exp = key
            if not value['combine_perturbations']:
                exp = "{}.Perturbation{}".format(key, perturbation)
            for sim in range(sims):
                for rep in range(reps):
                    experiments.append('{}Sim{}Rep{}'.format(exp, sim + 1,
                                                             rep + 1))
    return experiments

def get_experiment_ids(run_ids):
    """
    Get the list of experiment ids from the list of simulation ids
    
    Parameters
    ----------
    run_ids : list
        A list of individual simulation ids. Output from get_simulation_ids().
    
    Returns
    -------
    list
        Experiment ids with Sim and Rep numbers ripped from the original
        simulation id.
    """
    ids = set()
    sim_reg = re.compile('Sim*[0-9]')
    rep_reg = re.compile('Rep*[0-9]')
    for each in run_ids:
        new_id = rep_reg.sub('', sim_reg.sub('', each))
        new_id = new_id.replace('.', '')
        ids.add(new_id)
    return list(ids)