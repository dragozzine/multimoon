import numpy as np
import mm_priors as prior
import pandas as pd
import mm_param

"""
Inputs:
1) fit_array, the array of all the fitted parameters

Outputs:
1) log_likelihood, the log likelihood of the parameters with the priors

"""
def log_likelihood(params, obsdf, runprops, fitarray_to_params_dict):
    # assuming Guassian independent observations log-likehood = -1/2 * chisqure

    # DS TODO: convert parameters array to parameters dataframe
    
    #We currently do not have a defined constraints dictionary anywhere, so currently we'll define it as empty
    constraints = {}
    
    paramdf = mm_param.from_fit_array_to_param_df(params, constraints, fitarray_to_params_dict)

    lh = mm_chisquare(paramdf,obsdf)*-0.5
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
1)The Parameters dataframe
2) The Observation Dataframe

Outputs:
1) The chi-squared number of the likelihood
"""
# calculates the chi-square for parameters given observations
def mm_chisquare(paramdf, obsdf):
    
    
    # SP TODO: fill this in
    # use parameters dataframe with one set of parameters and observation times to call SPINNY to get the model dataframe
    # SPINNY returns (full?) state vector of each object (in primaricentric ecliptic J2000 coordinates) at each time
    
    # Are there "holes" in the dataframe when no value is initialized?
    # what are the dimensions of this dataframe? Is it only one set of parameters?
    
    time_arr = obsdf['time'].values.flatten() # gets an array of observation times from the obs dataframe
    T = len(time_arr)
    
    vec_df = spinny_vector(paramdf, time_arr) # returns a dataframe of state vectors for each body in the system
    
    
    for t in range(0,T):
        # what happens if the sun isn't included in the model? Can you still get heliocentric position?
        # which numbers do I want from the obsdf?
        
        convert_ecl_rel_pos_to_geo_rel_ast(ecl_rel_pos, obj_rel_pos, rel_moon):
    
    # use mm_relast code to turn positions into geocentric relative astrometry

    # calculate chisquare = sum [ (model-obs)/err ] ^2 

    # return chisquare
    
    return 1
