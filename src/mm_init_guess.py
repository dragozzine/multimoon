import numpy as np
import pandas as pd
"""
This file will have the user input their initial guess for where to start in the emcee program
"""

"""
Function to generate an initial guess for the emcee function
Inputs: a filename or datafrme in a parameter-generator format, int nwalkers
Output: "init_guess" dataframe in parameter format, nwalkers rows
"""
def mm_init_guess(runprops, nwalkers):
    """This function will produce the initial guess used in nPSF.
    
    Input: 
    
    runprops- All run properties for the code. Will include
        the name of the init_guess datafraem csv file.
    nwalkers - number of walkers. Added by Rochelle to make
        compatible with run_props in npsf_run.py 4/17/20
    
    Returns:
    
    params_df - A parameters dataframe with the same column names
        as start_guess_df and nwalker rows drawn from the 
        distribution.
    """
    filename = "../data/" + runprops("objname") + "/" + objname +"init_guess.csv"
    start_guess_df = pd.read_csv(filename)
    arrSet = start_guess_df.as_matrix()
      
    # Some code to help us get the names for the columns.
    name_dict = {0:"Junk"}
    n = 0
    dist_arr = []
    for col in start_guess_df.columns:
        name_dict[n] = col
        infos = start_guess_df[col].as_matrix()
        mean1, stdev1 = infos[0],infos[1]
        if n == 0:
            dist_arr = np.random.normal(mean1,stdev1,nwalkers)
        else:
            dist_arr = np.vstack((dist_arr,np.random.normal(mean1,stdev1,nwalkers)))
        n += 1
    
    dict_vals = []
    for x in name_dict:
        dict_vals.append(name_dict[x])
        
    dist_arr = np.transpose(dist_arr)
    
    params_df = pd.DataFrame(dist_arr, columns = [dict_vals])
    
    return params_df
    
    
"""
Function to generate an initial guess for the emcee function
Inputs: a filename or datafrme in a parameter-generator format, int nwalkers
Output: "init_guess" dataframe in parameter format, nwalkers rows
"""     
def get_init_guess_from_multirow_startguess(dataframe, nwalkers):
    if len(startguess) > nwalkers:
        
    if len(startguess) < nwalkers:
        
    return 1
        
        
def get_init_guess_from_singlerow_startguess_rough():
    return 1
    
def get_init_guess_from_prior():
    return 1
    