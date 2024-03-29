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

    # Reading in/getting necessary data
    filename = runprops.get("init_filename")
    nwalkers = runprops.get("nwalkers")
    fix_float = runprops.get("float_dict")

    start_guess_df = pd.read_csv(filename,sep=',',index_col=0)
    start_guess_df = start_guess_df.transpose()
    arrSet = start_guess_df.to_numpy()

    # Some code to help us get the names for the columns.
    name_dict = {0:"Junk"}
    n = 0
    dist_arr = []
    
    cut_df = pd.DataFrame()
    
    for i in range(runprops.get('numobjects')):
        for col in start_guess_df.columns:
            if '_'+str(i+1) in col or 'offset' in col:
                cut_df[col] = start_guess_df[col]
        

    for col in start_guess_df.columns:
        # In the future this should be moved to spinny stuff???
        if 'period' in col:
            for i in range(runprops.get('numobjects')):
                if str(i+1) in col:
                    name_dict[n] = 'sprate_'+str(i+1)
            name_dict[n] = col
            infos = 2*np.pi/start_guess_df[col].to_numpy()/60/60
            mean1, stdev1 = infos[0],infos[1]
        else:
            name_dict[n] = col
            infos = start_guess_df[col].to_numpy()
            mean1, stdev1 = infos[0],infos[1]
        
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
