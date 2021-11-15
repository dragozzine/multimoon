import pandas as pd
import numpy as np
import sys
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

    total_df_names = dataframe.columns
    
    fit_names = []

    for i in range(0,runprops.get('numobjects')):
        if runprops.get('lockspinanglesflag') == True:
            #print(runprops.get('dynamicstoincludeflags'))
            if int(runprops.get('dynamicstoincludeflags')[i]) != 0:
                #print(dataframe[['spaop_'+str(i+1)]].values)
                #print(dataframe[['aop_'+str(i+1)]].values)
                if int(runprops.get('dynamicstoincludeflags')[i]) == 2:
                    dataframe[['spaop_'+str(i+1)]] = dataframe[['aop_2']].values
                dataframe[['spinc_'+str(i+1)]] = dataframe[['inc_2']].values
                dataframe[['splan_'+str(i+1)]] = dataframe[['lan_2']].values
                if fix_float_dict.get('spaop_'+str(i+1)) == 1:
                    print('Since you have chosen to lock the spin angles, please change the spaop_'+str(i+1)+' variable in the float_dict to be fixed.')
                    sys.exit()
                if fix_float_dict.get('spinc_'+str(i+1)) == 1:
                    print('Since you have chosen to lock the spin angles, please change the spinc_'+str(i+1)+' variable in the float_dict to be fixed.')
                    sys.exit()
                if fix_float_dict.get('splan_'+str(i+1)) == 1:
                    print('Since you have chosen to lock the spin angles, please change the splan_'+str(i+1)+' variable in the float_dict to be fixed.')
                    sys.exit()
                #'''
    
    if runprops.get('transform'):
        if runprops.get('numobjects') > 3:
            print('Warning: Only masses 1-3 will be used in the transformations for now. Future work can be done later to increase this')
        if fix_float_dict.get('mass_1') == 1 and fix_float_dict.get('mass_2') == 1:
            if fix_float_dict.get('mass_3') == 1:
                dataframe[['mass_2']] = np.array(dataframe[['mass_1']])+np.array(dataframe[['mass_2']])
                dataframe[['mass_3']] = np.array(dataframe[['mass_3']])+np.array(dataframe[['mass_2']])
                
                fit_names.append('mass1+2')
                fit_names.append('mass1+2+3')
            else:
                dataframe[['mass_2']] = np.array(dataframe[['mass_1']])+np.array(dataframe[['mass_2']])
                fit_names.append('mass1+2')
        
        for i in range(runprops.get('numobjects')-1):
            pomega = np.array(dataframe[['aop_'+str(i+2)]])+np.array(dataframe[['lan_'+str(i+2)]])
            Lambda = pomega + np.array(dataframe[['mea_'+str(i+2)]])
            
            
            fit_names.append('lambda_'+str(i+2))
            fit_names.append('pomega_'+str(i+2))
            
            if fix_float_dict.get('lan_'+str(i+2)) == 1 and fix_float_dict.get('aop_'+str(i+2)) == 1:
                dataframe[['aop_'+str(i+2)]] = pomega
            
            if fix_float_dict.get('mea_'+str(i+2)) == 1 and fix_float_dict.get('aop_'+str(i+2)) == 1:                   
                dataframe[['mea_'+str(i+2)]] = Lambda
                        
            if fix_float_dict.get('ecc_'+str(i+2)) == 1 and fix_float_dict.get('aop_'+str(i+2)) == 1:
                ecc = np.array(dataframe[['ecc_'+str(i+2)]])
                pomega_rad = np.array(dataframe[['aop_'+str(i+2)]])*np.pi/180
            
                dataframe[['ecc_'+str(i+2)]] = np.array(ecc)*np.sin(pomega_rad)
                dataframe[['aop_'+str(i+2)]] = np.array(ecc)*np.cos(pomega_rad)
               
            if fix_float_dict.get('inc_'+str(i+2)) == 1 and fix_float_dict.get('lan_'+str(i+2)) == 1:
                inc = np.array(dataframe[['inc_'+str(i+2)]])*np.pi/180
                lan = np.array(dataframe[['lan_'+str(i+2)]])*np.pi/180
                 
                dataframe[['inc_'+str(i+2)]] = np.tan(inc/2)*np.sin(lan)
                dataframe[['lan_'+str(i+2)]] = np.tan(inc/2)*np.cos(lan)
        
        for i in range(runprops.get('numobjects')):
            if fix_float_dict.get('spinc_'+str(i+1)) == 1 and fix_float_dict.get('splan_'+str(i+1)) == 1:
                spinc = np.array(dataframe[['spinc_'+str(i+1)]])*np.pi/180
                splan = np.array(dataframe[['splan_'+str(i+1)]])*np.pi/180
                print('spinc ',i,' ', spinc)
                print('splan ',i,' ', splan)
                a = np.cos(spinc/2)*np.sin(splan)
                b = np.cos(spinc/2)*np.cos(splan)
                dataframe[['spinc_'+str(i+1)]] = a
                dataframe[['splan_'+str(i+1)]] = b
                #print('a ', a)
                #print('b ', b)
                
                #if (a[:]>np.sin(splan[:])).any():
                #    print("There is a greater a than splan")
                #    print("spinc:", a)
                #    print("splan:", splan)
                
                #dataframe[['spinc_'+str(i+1)]] = np.tan(spinc/2)*np.sin(splan)
                #dataframe[['splan_'+str(i+1)]] = np.tan(spinc/2)*np.cos(splan)
    
    num = 0
    fit_scale = dataframe.iloc[0]
    #print(fit_scale)
    fit_scale = fit_scale.to_frame().transpose()
    #Scale every column down by the values in the first row.
    
    for col in dataframe.columns:
        if fit_scale[col][0] != 0.0:
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

    for col in fit_scale.columns:
        fit_scale.rename(columns={col: col[0]}, inplace=True)
    j = 1
    for i in runprops.get('dynamicstoincludeflags'):
        if int(i) > 0:
            fit_names.append('period_'+str(j))
        j = j+1
    if int(runprops.get('dynamicstoincludeflags')[0]) > 0:
        for i in range(runprops.get('numobjects')-1):
            fit_names.append('sat_spin_inc_'+str(i+2))
    #print(fit_names)
    #fit_names = np.array(fit_names)
    return float_arr, float_names, fixed_df, total_df_names, fit_scale, fit_names
    
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
    #if runprops.get('includesun') == 1:
    #    np.delete(float_array, np.s_[0:5:1],0)
    #    np.delete(float_names, np.s_[0:5:1],0)
        
    
    Index = range(len(fixed_df.index))
    float_df = pd.DataFrame(data = [float_array],index = Index, columns = float_names)
    param_df = pd.DataFrame()
    fit_params = pd.DataFrame()
    
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
        undo_pomega = np.zeros(runprops.get('numobjects')-1)
        undo_pomega[:] = False
        undo_masses = np.zeros(runprops.get('numobjects')-1)
        undo_masses[:] = False
        undo_spin = np.zeros(runprops.get('numobjects'))
        undo_spin[:] = False
        
        for i in range(runprops.get('numobjects')-1):
            if 'ecc_'+str(i+2) in float_names and 'aop_'+str(i+2) in float_names:
                undo_ecc_aop[i] = True
                        
            if 'inc_'+str(i+2) in float_names and 'lan_'+str(i+2) in float_names:
                undo_inc_lan[i] = True
                
            if 'mea_'+str(i+2) in float_names and 'aop_'+str(i+2) in float_names:
                undo_lambda[i] = True
                
            if 'lan_'+str(i+2) in float_names and 'aop_'+str(i+2) in float_names:
                undo_pomega[i] = True
        
        for i in range(runprops.get('numobjects')):
            if 'spinc_'+str(i+1) in float_names and 'splan_'+str(i+1) in float_names:
                undo_spin[i] = True
        
        if 'mass_1' in float_names and 'mass_2' in float_names:
            if 'mass_3' in float_names and runprops.get('numobjects') > 2:
                undo_masses[1] = True
            else:
                undo_masses[0] = True
        
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
        #print(fit_scale.columns)
        for col in fit_scale.columns:
            param_col = col
            if type(col) != str:
                param_col = col[0]
            param_df[param_col] = param_df[param_col]*fit_scale.get(col)
            
    param_df = param_df.iloc[[0]]
    #print(param_df)
    #print(fit_scale)
    
    
    if runprops.get('transform'):
        
        if runprops.get('numobjects') > 2 and undo_masses[1]:
            mass_3 = np.array(param_df['mass_3'])
            mass_2 = np.array(param_df['mass_2'])
            mass_1 = np.array(param_df['mass_1'])
            
            fit_params['mass1+2+3'] = mass_3
            fit_params['mass1+2'] = mass_2
            
            param_df['mass_3'] = mass_3-mass_2
            param_df['mass_2'] = mass_2-mass_1   
        elif undo_masses[0]:
            mass_2 = np.array(param_df['mass_2'])
            mass_1 = np.array(param_df['mass_1'])
            
            fit_params['mass1+2'] = mass_2
            param_df['mass_2'] = mass_2-mass_1
        
        for i in range(runprops.get('numobjects')-1):
            if undo_ecc_aop[i]:
                a = np.array(param_df['ecc_'+str(i+2)])
                b = np.array(param_df['aop_'+str(i+2)])
                
                pomega = np.arctan2(a,b)*180/np.pi
                                
                if pomega < 0:
                    pomega = pomega%360
                param_df['aop_'+str(i+2)] = pomega
                param_df['ecc_'+str(i+2)] = a/np.sin(pomega*np.pi/180)
            
            if undo_inc_lan[i]:
            
                a = np.array(param_df['inc_'+str(i+2)])
                b = np.array(param_df['lan_'+str(i+2)])
                
                lan = np.arctan2(a,b)*180/np.pi
                if lan < 0:
                    lan = lan%360
                
                c = np.sin(lan*np.pi/180)
                
                inc = np.arctan2(a,c)*2*180/np.pi
                
                if inc < 0:
                    inc = inc%180

                param_df['inc_'+str(i+2)] = inc
                param_df['lan_'+str(i+2)] = lan
                           
            if undo_lambda[i]:
                fit_params['lambda_'+str(i+2)] = np.array(param_df['mea_'+str(i+2)])
                mea = np.array(param_df['mea_'+str(i+2)])-np.array(param_df['aop_'+str(i+2)])
                if mea < 0:
                    mea = mea%360
                elif mea > 360:
                    mea = mea%360
                param_df['mea_'+str(i+2)] = mea
                
            
            if undo_pomega[i]:
                fit_params['pomega_'+str(i+2)] = np.array(param_df['aop_'+str(i+2)])
                aop  = np.array(param_df['aop_'+str(i+2)])-np.array(param_df['lan_'+str(i+2)])
                if aop < 0:
                    aop = aop%360
                elif aop > 360:
                    aop = aop%360
                param_df['aop_'+str(i+2)] = aop
                
        for i in range(runprops.get('numobjects')):      
            if undo_spin[i]:
            
                a = np.array(param_df['spinc_'+str(i+1)])
                b = np.array(param_df['splan_'+str(i+1)])
                
                splan = np.arctan2(a,b)*180/np.pi
                #splan = np.arctan(a/b)*180/np.pi
                if splan < 0:
                    splan = splan%360
                
                c = np.sin(splan*np.pi/180)
                
                if a/c > 1 or a/c < -1:
                    print('a/c is', a/c, ', causing the error')
                    print('a ',a)
                    print('b ',b)
                    print('c ',c)
                    print('splan ', splan)
                    param_df['spinc_'+str(i+1)] = -1
                    param_df['splan_'+str(i+1)] = -1

                else:
                    #print(np.arccos(a/c))
                    spinc = np.arccos(a/c)*2*180/np.pi
                    #spinc = np.arctan2(a,c)*2*180/np.pi
                
                    if spinc < 0:
                        spinc = spinc%180
                    
                    param_df['spinc_'+str(i+1)] = spinc
                    param_df['splan_'+str(i+1)] = splan
                    #print('spinc_new ',i,' ', spinc)
                    #print('splan_new ',i,' ', splan)
                
            if int(runprops.get('dynamicstoincludeflags')[0]) > 0:
                spinc1=np.deg2rad(np.array(param_df['spinc_1']))
                splan1=np.deg2rad(np.array(param_df['splan_1']))
                for i in range(runprops.get('numobjects')-1):
                    inc = np.deg2rad(np.array(param_df['inc_'+str(i+2)]))
                    lan = np.deg2rad(np.array(param_df['lan_'+str(i+2)]))
                    mutualinc = np.arccos( np.cos(spinc1)*np.cos(inc) + np.sin(spinc1)*np.sin(inc)*np.cos(splan1 - lan) )
                    mutualinc = np.rad2deg(mutualinc)
                    fit_params['spin_sat_inc_'+str(i+2)] = mutualinc
        
        N = runprops.get('numobjects')
        '''
        if undo_spin[N-1]:
            a = np.array(param_df['spinc_'+str(N)])
            b = np.array(param_df['splan_'+str(N)])
                
            splan = np.arctan2(a,b)*180/np.pi
            if splan < 0:
                splan = splan%360
                
            c = np.sin(splan*np.pi/180)
                
            if a/c > 1 or a/c < -1:
                print('a/c is', a/c, ', causing the error')
                param_df['spinc_'+str(i+1)] = -1
                param_df['splan_'+str(i+1)] = -1
            else:
                #print(np.arccos(a/c))
                spinc = np.arccos(a/c)*2*180/np.pi
                #spinc = np.arctan2(a,c)*2*180/np.pi
            
                if spinc < 0:
                    spinc = spinc%180

                param_df['spinc_'+str(i+1)] = spinc
                param_df['splan_'+str(i+1)] = splan
        '''        
    if runprops.get('lockspinanglesflag') == True:
        param_df['spinc_1'] = param_df['inc_2']
        param_df['splan_1'] = param_df['lan_2']

    return param_df, fit_params

