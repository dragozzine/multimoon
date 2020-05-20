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
    return 1
    
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