import json

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
            exp_pert = "{}.Perturbation{}".format(key, perturbation)
            for sim in range(sims):
                for rep in range(reps):
                    experiments.append('{}Sim{}Rep{}'.format(exp_pert, sim +1,
                                                            rep + 1))
    return experiments