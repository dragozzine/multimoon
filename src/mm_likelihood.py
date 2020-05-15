"""
Inputs:
1) fit_array, the array of all the fitted parameters

Outputs:
1) log_likelihood, the log likelihood of the parameters with the priors

"""
def log_likelihood(obsDF):
    
    return 1

"""
Inputs:
1)The Model Dataframe
2) The Observation Dataframe

Outputs:
1) The chi-squared number of the likelihood
"""
def mm_chisquare(modeldata, obsdata):
    
    return 1


def log_probability(fit_array, args, kwargs):
    import mm_prior
    
    mm_prior.mm_prior(priors,param)
    
    
    llhood = log_likelihood(obsDF, paramDF)
    
    