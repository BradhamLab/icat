import inspect
import numpy as np
import collections

def check_kws(reference_dict, new_dict, name):
    if not isinstance(new_dict, dict):
        raise ValueError("Expected dictionary of keyword arguments for "
                         "`{}`. Received {}.".format(name, type(new_dict)))
    for key, item in new_dict.items():
        if key not in reference_dict.keys():
            raise ValueError("Unsupported keyword argument `{}` for "
                             "{} keywords.".format(key, name))
        new_dict[key] = item
    return new_dict


def check_matching_genes(ref, new):
    return set(ref.var.index.values).difference(new.var.index.values) == 0


def get_default_kwargs(func, ignore_params=[]):
    params = inspect.signature(func).parameters
    kwargs = {x:params[x].default for x in params if x not in ignore_params}
    return kwargs


def check_np_castable(obj, name):
    """Check whether an object is castable to a numpy.ndarray."""
    if not isinstance(obj, np.ndarray):
        try:
            obj = np.array(obj)
        except:
            raise ValueError("Expected numpy.ndarray castable object for "
                            "`{}`. Got {}.".format(obj, type(obj)))
    return obj


def __evaluate_key(key, sep):
    if not isinstance(key, str):
        raise ValueError("Keys must be strings for easy concatentation.")
    if sep in key:
        raise ValueError("Cannot have `{}` in dictionary keys.".format(sep))


def flatten_dict(d, parent_key='', sep='.'):
    """
    Flatten a dictionary containing nested dictionaries.
    
    Parameters
    ----------
    d : dict
        Dictionary to flatten
    parent_key : str, optional
        Key in parent dictionary pointing to `d`. The default is '', which
        assumes `d` is the highest level nested dictionary.
    sep : str, optional
        String value to separate child and parent keys. The default is '.',
        which will place a '.' between each key. All parent and child keys
        will be assessed to ensure they do not contain a `sep` character;
        therefore, `sep` should be set to a delimiter not present in current
        keys.
    
    Returns
    -------
    dict
        Flattened dictionary with parent and child keys separted by `sep`.

    References
    ----------

    Taken shamelessly from here:
        https://stackoverflow.com/questions/6027558/flatten-nested-python-dictionaries-compressing-keys
    """

    items = []
    for k, v in d.items():
        __evaluate_key(k, sep)
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)