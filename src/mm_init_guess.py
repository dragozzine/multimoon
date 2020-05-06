import numpy as np

def mm_init_guess(runprops, walkers):
    ndim = 2
    p0 = np.random.random((ndim, walkers))
    return p0
"""
This file will have the user input their initial guess for where to start in the emcee program
"""

"""
Function to generate an initial guess for the emcee function
Inputs: a filename or dataframe in a parameter-generator format, int nwalkers
Output: "init_guess" dataframe in parameter format, nwalkers rows
"""
def gen_init_guess_from_dist(dataframe, nwalkers):
    
    return 1
    
    
"""
Function to generate an initial guess for the emcee function
Inputs: a filename or datafrme in a parameter-generator format, int nwalkers
Output: "init_guess" dataframe in parameter format, nwalkers rows
"""     
def get_init_guess_from_multirow_startguess(dataframe, nwalkers):
    c = 0
    if len(startguess) > nwalkers:
        c=1
    if len(startguess) < nwalkers:
        c=1
    return 1
        
        
def get_init_guess_from_singlerow_startguess_rough():
    return 1
    
def get_init_guess_from_prior():
    return 1
    