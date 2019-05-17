import inspect

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