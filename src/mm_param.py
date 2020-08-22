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
    
    for i in range(runprops.get('numobjects')):
        if runprops.get('lockspinanglesflag') and runprops.get('dynamicstoincludeflags')[i+1] != 0:
            dataframe['spaop_'+str(i+1)] = dataframe['aop_'+str(i+1)]
            dataframe['spinc_'+str(i+1)] = dataframe['inc_'+str(i+1)]
            dataframe['splan_'+str(i+1)] = dataframe['lan_'+str(i+1)]
            if fix_float_dict.get('spaop_'+str(i+1)) == 1:
                print('Since you have chosen to lock the spin angles, please change the spaop_'+str(i+1)+' variable in the float_dict to be fixed.')
                sys.exit()
            if fix_float_dict.get('spinc_'+str(i+1)) == 1:
                print('Since you have chosen to lock the spin angles, please change the spinc_'+str(i+1)+' variable in the float_dict to be fixed.')
                sys.exit()
            if fix_float_dict.get('splan_'+str(i+1)) == 1:
                print('Since you have chosen to lock the spin angles, please change the splan_'+str(i+1)+' variable in the float_dict to be fixed.')
                sys.exit()
    
    num = 0
    #Scale every column down by the values in the first row.
    for col in dataframe.columns:
        dataframe[col] = dataframe[col]/fit_scale[col][0]
        num = num+1
    
    key_list = list(fix_float_dict.keys()) 
    val_list = list(fix_float_dict.values())
    
    fixed_df = pd.DataFrame(index = range(len(dataframe.index)))
    float_df = pd.DataFrame()
    float_names = []
    num = 0
    float_array = np.array([])
    if len(key_list) == 0:
        float_array = dataframe.to_numpy()
    else:
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
        
        for i in range(runprops.get('numobjects')-1):
            
            if fix_float_dict.get('mass_'+str(i+2)) == 1 and fix_float_dict.get('aop_'+str(i+2)) == 1:
                float_df['mass_'+str(i+2)] = float_df['mass_'+str(i+2)]+float_df['aop_'+str(i+2)]
            
            if fix_float_dict.get('ecc_'+str(i+2)) == 1 and fix_float_dict.get('aop_'+str(i+2)) == 1:
                ecc = float_df['ecc_'+str(i+2)]
                aop = float_df['aop_'+str(i+2)]
            
                float_df['ecc_'+str(i+2)] = ecc*np.cos(aop)
                float_df['aop_'+str(i+2)] = ecc*np.sin(aop)
    
            if fix_float_dict.get('inc_'+str(i+2)) == 1 and fix_float_dict.get('lan_'+str(i+2)) == 1:
                inc = float_df['inc_'+str(i+2)]
                lan = float_df['lan_'+str(i+2)]
            
                float_df['inc_'+str(i+2)] = np.tan(inc/2)*np.sin(lan)
                float_df['lan_'+str(i+2)] = np.tan(inc/2)*np.cos(lan)
                          
    
    
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
def from_fit_array_to_param_df(float_array, float_names, fixed_df, total_df_names, fit_scale, names_dict, runprops):
    
    #First, turn the float_array back into  dataframe with the column names given
    Index = range(len(fixed_df.index))
    float_df = pd.DataFrame(data = [float_array],index = Index, columns = float_names)
    
    param_df = pd.DataFrame()
    
    
    if len(fixed_df) == 0:
        param_df = float_df
    else:
    #Recombine the float and fixed dataframes
        undo_ecc_aop = np.zeros(runprops.get('numobjects')-1)
        undo_ecc_aop[:] = False
        undo_inc_lan = np.zeros(runprops.get('numobjects')-1)
        undo_inc_lan[:] = False
        undo_lambda = np.zeros(runprops.get('numobjects')-1)
        undo_lambda[:] = False
        
        for i in range(runprops.get('numobjects')-1):
            if 'ecc_'+str(i+2) in float_names and 'aop_'+str(i+2) in float_names:
                undo_ecc_aop[i] = True
            
            if 'inc_'+str(i+2) in float_names and 'lan_'+str(i+2) in float_names:
                undo_inc_lan[i] = True
                
            if 'mass_'+str(i+2) in float_names and 'aop_'+str(i+2) in float_names:
                undo_lambda[i] = True
                
                
        for i in total_df_names:
            name = i[0]
            if name in fixed_df:
                value = fixed_df[name].values.tolist()
                param_df[name] = value
            elif name in float_df:
                value = float_df[name].values.tolist()
                param_df[name] = value
        
        names_df = pd.DataFrame.from_dict(names_dict,orient='index')
        names_df = names_df.transpose()
        
        for col in names_df.columns:
            param_df[col] = names_df[col][0]
    
      
        #Now unfit all of the variables by multipliyng each column by its fit variable.
        #print(param_df)
        for col in fit_scale.columns:
            param_df[col[0]] = param_df[col[0]]*fit_scale[col][0]
            
    param_df = param_df.iloc[[0]]
    
    for i in range(runprops.get('numobjects')-1):
        if undo_ecc_aop[i]:
            param_df['aop_'+str(i+2)] = np.arctan(np.array(param_df['aop_'+str(i+2)])/np.array(param_df['ecc_'+str(i+2)]))
            param_df['ecc_'+str(i+2)] = param_df['ecc_'+str(i+2)]/np.sin(np.array(param_df['aop_'+str(i+2)]))
            
        if undo_inc_lan[i]:
            
            sinlan = param_df['inc_'+str(i+2)]
            coslan = param_df['lan_'+str(i+2)]
            
            inc = np.arctan(np.array(sinlan)/np.array(coslan))
            param_df['inc_'+str(i+2)] = inc
            
            param_df['lan_'+str(i+2)] = np.arcsin(np.array(sinlan)/np.tan(inc/2))
            
        if undo_lambda[i]:
            param_df['mass_'+str(i+2)] = param_df['mass_'+str(i+2)]-param_df['aop_'+str(i+2)]
            
                              
           
    return param_df
