"""
Function to convert the parameter dataframe to a scaled and fitted array.
Inputs: 
1) The Parameters dataframe
2) The fix/float/constrain constraints dictionary
3) The dictionary describing the scale of each element

Outputs:
1) The fitted array of parameters
2) The dictionary of the param fit data
"""
def from_param_df_to_fit_array(dataframe, constraints={}, param_to_fit_scale={}):
    
    cols = dataframe.columns
    numWalkers = len(dataframe.index)
    
    param_arr = []
    
    for i in range(numWalkers):
        num = 0
        for j in cols:
            param_arr[num][i] = dataframe[j][i]
            num = num+1
    
    return param_arr
    
"""
Function to convert a fitted array into the parameter dataframe
Inputs: 
1) The fitted array
2) The fix/float/constrain constraints dictionary
3) The dictionary describing the scale of each element

Outputs:
1) Dataframe in parameter format
"""
def from_fit_array_to_param_df(fit_array, constraints={}, param_to_fit_scale={}):
    import runprops
    
    dict_vals = runprops.runprops.get("float_dict").keys()
  
    # We make the dataframe
    param_df = pd.DataFrame(fit_array,
                            columns = [dict_vals])
    
    return param_df