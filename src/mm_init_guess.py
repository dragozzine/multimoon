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
def mm_init_guess(runprops):
    """This function will produce the initial guess used in multimoon.
    
    Input: 
    
    runprops- All run properties for the code. Will include
        the name of the init_guess dataframe csv file.
    
    Returns:    
    params_df - A parameters dataframe with the same column names
        as start_guess_df and nwalker rows drawn from the 
        distribution.
    """
    filename = runprops.get("init_filename")

    start_guess_df = pd.read_csv(filename,sep=',',index_col=0)

    start_guess_df = start_guess_df.transpose()

    arrSet = start_guess_df.to_numpy()

    nwalkers = runprops.get("nwalkers")
    fix_float = runprops.get("float_dict")
    
      
    # Some code to help us get the names for the columns.
    name_dict = {0:"Junk"}
    n = 0
    dist_arr = []
    
    cut_df = pd.DataFrame()
    
    for i in range(runprops.get('numobjects')):
        for col in start_guess_df.columns:
            if '_'+str(i+1) in col:
                cut_df[col] = start_guess_df[col]
        

    for col in cut_df.columns:
        name_dict[n] = col
        infos = cut_df[col].to_numpy()
        mean1, stdev1 = infos[0],infos[1]
        #print(col)
        #print(fix_float.get(col))
        
        if n == 0:
            if fix_float.get(col) == 0:
                dist_arr = mean1*np.ones(nwalkers)
            else:
                dist_arr = np.random.normal(mean1,stdev1,nwalkers)
        else:
            if fix_float.get(col) == 0:
                dist_arr = np.vstack((dist_arr, mean1*np.ones(nwalkers)))
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
         dum=0
    if len(startguess) < nwalkers:
         dum=0
    return 1
        
        
def get_init_guess_from_singlerow_startguess_rough():
    return 1
    
def get_init_guess_from_prior():
    return 1
    
