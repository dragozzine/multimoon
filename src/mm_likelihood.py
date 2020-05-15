import numpy as np
import mm_priors as prior
import pandas as pd

"""
Inputs:
1) fit_array, the array of all the fitted parameters

Outputs:
1) log_likelihood, the log likelihood of the parameters with the priors

"""
def log_likelihood(params, obsdf, runprops, fitarray_to_params_dict):
    lh = mm_chisquare(params,obsdf)*-0.5
    return lh


"""
Inputs:
1) params
2) runprops
3) fitarray_dict
4) 

Outputs:
1) log_probability, the log_likelihood plus the priors, which is the total probability

"""
def log_probability(params, runprops, fitarray_to_params_dict, obsdf):
    
    objname = runprops.get("objname")
    priors = pd.read_csv("../data/" +objname + "/" + objname + "_priors_df.csv", sep='\t',index_col=0)
    priors = priors.transpose()
    
    lp = prior.mm_priors(priors,params)
    
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, obsdf, runprops, fitarray_to_params_dict)


"""
Inputs:
1)The Model Dataframe
2) The Observation Dataframe

Outputs:
1) The chi-squared number of the likelihood
"""
def mm_chisquare(modeldata, obsdata):
    
    
    return 1
    