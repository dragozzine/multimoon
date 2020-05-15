import numpy as np
import mm_priors as prior
"""
Inputs:
1) fit_array, the array of all the fitted parameters

Outputs:
1) log_likelihood, the log likelihood of the parameters with the priors

"""
def log_likelihood(params, obsdf):
    lh = mm_chiquare(params,obsdf)*-0.5
    return lh
    


def log_probability(params, runprops, fitarray_to_params_dict, obsdf):
    lp = prior.mm_prior()
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, obsdf)


"""
Inputs:
1)The Model Dataframe
2) The Observation Dataframe

Outputs:
1) The chi-squared number of the likelihood
"""
def mm_chisquare(modeldata, obsdata):
    
    
    return 1
    