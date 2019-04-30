
class icat():
    """
    Model to Identify Clusters Across Treatments. 
    """

    def __init__(self, method='ncfs', method_kws=None):
        self.method = method
        self.method_kws = method_kws
    
    @property
    def method(self):
        return self._method
    
    @method.setter
    def method(self, value):
        if not isinstance(value, str):
            raise ValueError("Expected string for `method` parameter." 
                             "Received {}".format(value))
        value = value.lower()
        if value not in ['ncfs', 'lda', 'qda']:
            raise ValueError("Unsupported method: {}".format(value))
        self._method = value

    @property
    def method_kws(self):
        return self._method_kws
    
    @method_kws.setter
    def method_kws(self, value):
        if not isinstance(value, dict) and value is not None:
            raise ValueError("Expected dictionary of keyword arguments for "
                             "`method_kws`. Received {}.".format(type(value)))
        default_kws = {'ncfs': {'alpha': 0.01, 'metric': 'euclidean', 'reg': 3,
                                'sigma': 2, 'eta': 10e-6},
                       'lda': {'solver': 'eigen', 'shrinkage': 'auto',
                               'priors': None, 'n_components': None,
                               'store_covariance': False, 'tol': 0.0001},
                       'qda': {'solver': 'eigen', 'shrinkage': 'auto',
                               'priors': None, 'n_components': None,
                               'store_covariance': False, 'tol': 0.0001}}
        kws = default_kws[self.method]
        for key, item in value.items():
            if key not in kws.keys():
                raise ValueError("Unexpected keyword argument `{}` for method "
                                 "{}.".format(key, self.method))
            kws[key] = item 
        self._method_kws = kws