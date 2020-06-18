import pandas as pd
import numpy as np
"""
Function to convert the parameter dataframe to a scaled and fitted array.
Inputs: 
1) The Parameters dataframe
2) The fix/float/constrain constraints dictionary is inside runprops

Outputs:
1) The fitted array of floating parameters
2) The array of the names of each array
3) The dataframe of fixed values
4) The list of all column names in order for recombination
5) The original fit scale which will be needed later during recombination
"""
def from_param_df_to_fit_array(dataframe, runprops):
    
    fix_float_dict = runprops.get("float_dict")
    fit_scale = dataframe.iloc[0]
    fit_scale = fit_scale.to_frame().transpose()
    total_df_names = dataframe.columns
    
    num = 0
    #Scale every column down by the values in the first row.
    for col in dataframe.columns:
        dataframe[col] = dataframe[col]/fit_scale[col][0]
        num = num+1
    
    key_list = list(fix_float_dict.keys()) 
    val_list = list(fix_float_dict.values())
    
    fixed_df = pd.DataFrame()
    float_df = pd.DataFrame()
    float_names = []
    num = 0
    #Split the fixed and floating values into seperate dataframes
    for col in dataframe.columns:
        #If the value is fixed
        name = col[0]
        if fix_float_dict.get(col[0]) == 0:
            fixed_df[name] = dataframe[col]
        #If the value is floating
        elif fix_float_dict.get(col[0]) == 1:
            float_df[name] = dataframe[col]
            float_names.append(name)
        num = num+1
    
    float_arr = float_df.to_numpy()
    
    for col in fit_scale.columns:
        fit_scale.rename(columns={col: col[0]}, inplace=True)
    
    return float_arr, float_names, fixed_df, total_df_names, fit_scale
    
"""
Function to convert a fitted array into the parameter dataframe
Inputs: 
1) The fitted float array
2) names of each column in the float_array
3) Fixed dataframe of values
4) All names of parameters in order
5) The scale of fit variables

Outputs:
1) Dataframe in parameter format
"""
def from_fit_array_to_param_df(float_array, float_names, fixed_df, total_df_names, fit_scale, names_dict):
    
    #First, turn the float_array back into  dataframe with the column names given
    Index = range(len(fixed_df.index))
    float_df = pd.DataFrame(data = [float_array],index = Index, columns = float_names)
    
    param_df = pd.DataFrame()
    #Recombine the float and fixed dataframes
    for i in total_df_names:
        name = i[0]
        if name in fixed_df:
            value = fixed_df[name].values
            param_df[name] = value
        elif name in float_df:
            value = float_df[name].values
            param_df[name] = value
    
    names_df = pd.DataFrame.from_dict(names_dict,orient='index')
    names_df = names_df.transpose()
    
    for col in names_df.columns:
        param_df[col] = names_df[col][0]
    
    #Now unfit all of the variables by multipliyng each column by its fit variable.
    for col in fit_scale.columns:
        param_df[col[0]] = param_df[col[0]]*fit_scale[col]

    return param_df