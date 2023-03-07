
import numpy as np
import mm_priors as prior
import pandas as pd
import mm_param
import sys
sys.path.insert(1, 'mm_SPINNY/')
from mm_SPINNY.spinny_vector import generate_vector
import random
import mm_relast
from csv import writer
import os
import time
from scipy.stats import chi2

"""
Inputs:
1) fit_array, the array of all the fitted parameters
Outputs:
1) log_likelihood, the log likelihood of the parameters with the priors
"""
def log_likelihood(params, obsdf, runprops, geo_obj_pos):
    # assuming Gaussian independent observations log-likelihood = -1/2 * chisquare
    
    # Start by finding the chi-square and calculating log likelihood
    lh,residuals = mm_chisquare(params,obsdf, runprops, geo_obj_pos)
    if not np.isfinite(lh):
        return -np.inf, residuals
    lh = lh*-0.5

    # Robust statistics. Not working in MultiMoon v1.0
    if runprops.get("robust_stats"):
        rows = obsdf.shape[0]
        numObj = runprops.get("numobjects")
        ll_robust = 0
        jitter = params["jitter"].iloc[0]
        p_outlier = params["pbad"].iloc[0]

        names_dict = runprops.get("names_dict")
        names=[0 for i in range(numObj)]
        for i in range(0,numObj):
            names[i] = names_dict.get("name_"+str(i+1))

        # Loop over the all the data point for each object. This could be vectorized in the future.
        for j in range(1,numObj):
            for i in range(rows):
                # Extract uncertainties from the data
                lon_err = obsdf["DeltaLong_"+names[j]+"_err"][i]
                lat_err = obsdf["DeltaLat_"+names[j]+"_err"][i]
                
                # Calculate variance for the background model
                combinedlon_err = np.sqrt(lon_err**2 + jitter**2)
                combinedlat_err = np.sqrt(lat_err**2 + jitter**2)
                
                # Extract O-C from the residuals
                omc_lon = (residuals[2*(j-1)][i] * obsdf["DeltaLong_"+names[j]+"_err"][i])
                omc_lat = (residuals[2*(j-1)+1][i] * obsdf["DeltaLat_"+names[j]+"_err"][i])
                
                # Calculate likelihood for foreground model
                ll_fg_lon = -0.5 * (omc_lon/lon_err)**2 - np.log(lon_err)
                ll_fg_lat = -0.5 * (omc_lat/lat_err)**2 - np.log(lat_err)
                ll_fg = ll_fg_lon + ll_fg_lat
                
                # Calculate likelihood for background model
                ll_bg_lon = -0.5 * (omc_lon/combinedlon_err)**2 - np.log(combinedlon_err)
                ll_bg_lat = -0.5 * (omc_lat/combinedlat_err)**2 - np.log(combinedlat_err)
                ll_bg = ll_bg_lon + ll_bg_lat
                
                # Calculate the combined likelihood, including the penlaization for data rejection
                arg1 = np.log(1 - p_outlier) + ll_fg
                arg2 = np.log(p_outlier) + ll_bg
                if np.isnan(arg2):
                    print("args", arg1, arg2)
                    print("params", jitter, p_outlier)
                
                # Use logaddexp for numerical stability, as suggested by Foreman-Mackey
                ll_i = np.logaddexp(arg1, arg2)

                # Add total likelihood for data point to the summed likelihood
                if not np.isnan(ll_i):
                    ll_robust += ll_i

        return ll_robust, residuals

    # Return log likelihood and residuals arrays
    return lh, residuals


"""
Inputs:
1) params
2) runprops
3) fitarray_dict
4) 
Outputs:
1) log_probability, the log_likelihood plus the priors, which is the total probability
"""
def log_probability(float_params, float_names, fixed_df, total_df_names, fit_scale, runprops, obsdf, geo_obj_pos, best_llhoods):
    
    # Getting values from runprops
    objname = runprops.get("objectname")
    priorFilename = runprops.get('priors_filename')

    # Ensure that we are in the correct folder    
    if 'runs' in os.getcwd() or 'results' in os.getcwd():
        os.chdir('../../../src')

    # Loading in the priors
    priors = pd.read_csv(priorFilename, sep='\t',index_col=0)
    priors = priors.transpose()
    
    # Getting the names of each object
    name_dict = runprops.get("names_dict")
    if runprops.get('includesun') == True:
        name_dict['name_0'] = 'Sun'

    # Transform fitting units/parameters to actual units/parameters
    params,fit_params = mm_param.from_fit_array_to_param_df(float_params, float_names, fixed_df, total_df_names, fit_scale, name_dict, runprops)

    # Evaluate the priors    
    lp = prior.mm_priors(priors,params,runprops)

    # Output some things
    if runprops.get('verbose'):
        print('LogPriors: ',lp)
    if not np.isfinite(lp):
        return -np.inf

    # Get the log likelihood and residuals and calculate the combined probability
    log_likeli, residuals = log_likelihood(params, obsdf, runprops, geo_obj_pos)
    if not np.isfinite(log_likeli):
        return -np.inf
    llhood = lp + log_likeli

    # Begin machinery to output to the best likelihod file
    the_file = runprops.get('results_folder') + '/best_likelihoods.csv'

    #You will notice this differes from the regular runs way to save data
    #Since we are using mpi, we need to continually retrieve the best_likelihoods csv
    #I did the math with some intense testing, and found this will only slow down
    #a 1000 step system by 1 minute, which typically takes 2 hours, so there is not much slow down.
    
    #The best_llhoods dictionary keeps track of the best likelihood achieved by each individual processor, whic operates individually
    #from the othe rprocessors. If a processor/waler achieved a lieklihood that is better than the processor has achieved before,
    #we enter this if statement. 
    if llhood > best_llhoods.get("best_llhood") and runprops.get("is_mcmc") and runprops.get("updatebestfitfile") :
        best_llhoods['best_llhood'] = llhood
        best_llhoods['best_params'] = params.to_dict()
        best_csv = pd.read_csv(the_file, index_col=None)        
        
        #Here we detemrine the current best likelihood overall. This is saved in the best_likelihoods.csv, which is independent, and can be
        # accessed by any processor. This way the processors can all access the overall best likelihood value without overlapping very 
        #often.
        if len(best_csv.index) < 1:
            curr_best = -np.inf
        else:
            curr_best = best_csv.iloc[-1,0]
            
        if type(curr_best) == str:
            curr_best = 0
        if llhood > curr_best:
            chi_sq = llhood/(-0.5)            
            reduced_chi_sq = chi_sq/best_llhoods.get('deg_freedom')
            p_val = 1 - chi2.cdf(chi_sq, best_llhoods.get('deg_freedom'))
            with open(the_file, 'a+', newline='') as write_obj:
                csv_writer = writer(write_obj, delimiter = ',')
                thelist = params.head(1).values.tolist()[0]
                thelist.insert(0, lp)
                thelist.insert(0, reduced_chi_sq)
                thelist.insert(0,chi_sq)
                thelist.insert(0, p_val)
                thelist.insert(0, best_llhoods.get('deg_freedom'))
                
                thelist.insert(0, llhood)
                for i in range(runprops.get('numobjects')):
                    thelist.pop()
                if runprops.get('includesun'):
                    thelist.pop()
                for i in range(runprops.get("numobjects")-1):
                    thelist.append(residuals[2*(i-1)])
                    thelist.append(residuals[2*(i-1)+1])
                csv_writer.writerow(thelist)
    return llhood


"""
Inputs:
1)The Parameters dataframe
2) The Observation Dataframe
Outputs:
1) The chi-squared number of the likelihood
"""
# calculates the chi-square for parameters given observations
def mm_chisquare(paramdf, obsdf, runprops, geo_obj_pos, gensynth = False):

    # Get things from runprops
    numObj = runprops.get("numobjects")
    verbose = runprops.get("verbose")
    pd.set_option('display.max_columns', None)
    names = []
    for i in range(1,numObj+1):
        names.append('name_'+str(i))
        if not 'name_'+str(i) in paramdf.columns:
            print('The parameter name_' + str(i)+ ' is not found in the parameter dataframe.')
            sys.exit()
        
    # Sort observations file by times
    obsdf = obsdf.sort_values(by=['time'])

    # Get the times which need to be found
    time_arr = obsdf['time'].values.flatten()

    # Setting times relative to the epoch
    epoch = runprops.get("epoch_SJD")
    time_arr = time_arr - epoch

    # Start a timer to time how long this step takes
    begin = time.time()

    # Run spinny simulation inside of try except
    try:
        time_arr_sec = time_arr*86400
        vec_df = generate_vector(paramdf, time_arr_sec, runprops)
    except Exception as e:
        print('There was an error thrown within spinny:\n', e)
        rows = obsdf.shape[0]
        # Output inf if spinny has an error
        return np.inf, np.ones(((numObj-1)*2, rows))*10000

    # Making sure we have the right names for the objects
    names_dict = runprops.get("names_dict")
    names=[0 for i in range(numObj)]
    for i in range(0,numObj):
        names[i] = names_dict.get("name_"+str(i+1))

    # End timer
    end = time.time()

    # Check if the output vectors are primaricentric
    name_1 = "X_Pos_"+names[0]
    if (vec_df[name_1][0] != 0.0):
        print("Not primaricentric like I thought!")
        rows = obsdf.shape[0]
        return np.inf, np.ones(((numObj-1)*2, rows))*10000
    
    # Set up arrays
    Model_DeltaLong = np.zeros((numObj-1,len(time_arr)))
    Model_DeltaLat = np.zeros((numObj-1,len(time_arr)))
    if runprops.get('includesun') == 1:
        #print(vec_df)
        vec_df = vec_df.drop(['X_Pos_Sun', 'Y_Pos_Sun', 'Z_Pos_Sun', 'X_Vel_Sun', 'Y_Vel_Sun', 'Z_Vel_Sun'], axis=1)

    positionData = np.zeros((numObj*3,len(time_arr)))

    # Load positions into positionData
    for i in range(0,numObj):
        positionData[3*i] = vec_df["X_Pos_"+names[i]]
        positionData[3*i+1] = vec_df["Y_Pos_"+names[i]]
        positionData[3*i+2] = vec_df["Z_Pos_"+names[i]]

    # Changing to geocentric positions
    obs_to_prim_pos = [positionData[0]+geo_obj_pos['x'].tolist(),positionData[1]+geo_obj_pos['y'].tolist(),positionData[2]+geo_obj_pos['z'].tolist()]

    # Converting x,y,z coordinates to relative astrometry
    for i in range(1,numObj):
        prim_to_sat_pos = [positionData[i*3],positionData[i*3+1],positionData[i*3+2]]
        Model_DeltaLong[i-1], Model_DeltaLat[i-1] = mm_relast.convert_ecl_rel_pos_to_geo_rel_ast(obs_to_prim_pos, prim_to_sat_pos)
    
    '''
    # We are commenting this out until MultiMoon 2.0. We will release this as a feature then.
    # Putting in photocenter-braycenter offset for hidden objects in >2 object systems
    if runprops.get('photo_offset'):
        
        mass_ratio = paramdf['mass_2'][0]/paramdf['mass_1'][0]
        
        f_val = paramdf['f_val_1'][0]
        
        bright_ratio = f_val*mass_ratio**(2/3)
        
        rel_pos_lat = Model_DeltaLat[0,:]
        rel_pos_long = Model_DeltaLong[0,:]
        

        delta_offset_lat = bright_ratio*rel_pos_lat
        delta_offset_long = bright_ratio*rel_pos_long
        
        Model_DeltaLat = Model_DeltaLat - delta_offset_lat
        Model_DeltaLong = Model_DeltaLong - delta_offset_long
    '''
    # Adding in center of mass center of light offsets    
    if runprops.get("com_offset"):
        Model_DeltaLong = Model_DeltaLong + paramdf["long_offset"].iloc[0]
        Model_DeltaLat = Model_DeltaLat + paramdf["lat_offset"].iloc[0]

    # Outputting the Model_DeltaLong and Lat if gensynth flag is included in function call
    if gensynth:
        print("Returning the Model_DeltaLong and Lat dataframes for use in synthetic astrometry.")
        return Model_DeltaLong, Model_DeltaLat, obsdf

    # Now we have model delta Long and delta Lat for each object and each time 
    rows = obsdf.shape[0]

    # Setting up storage arrays
    residuals = np.zeros(((numObj-1)*2, rows))
    get_residuals = runprops.get("get_resid")
    delta_offset = 0
        
    # Calculating the residuals
    for i in range(rows):
        for j in range(1,numObj):

            residuals[2*(j-1)][i] = ((Model_DeltaLong[j-1][i]-obsdf["DeltaLong_"+names[j]][i])/obsdf["DeltaLong_"+names[j]+"_err"][i])
            residuals[2*(j-1)+1][i] = ((Model_DeltaLat[j-1][i]-obsdf["DeltaLat_"+names[j]][i])/obsdf["DeltaLat_"+names[j]+"_err"][i])
           
            if verbose:
                print("i,j,model,obs,err")
                print(i, j, Model_DeltaLong[j-1][i], obsdf["DeltaLong_"+names[j]][i], obsdf["DeltaLong_"+names[j]+"_err"][i])

                                      
    # Calculating the chisqaure for each observation
    chisquares = residuals**2
    
    # Calculating chi-square for all objects at all times
    chisq_tot = np.zeros(2*numObj)
    for i in range(0,2*numObj-2):
        chisq_tot[i]=np.nansum(chisquares[i])
        
    chisquare_total = np.nansum(chisq_tot)

    if verbose:
        print("chisq_tot, chisquare_total, residuals")
        print(chisq_tot, chisquare_total, residuals)

    return chisquare_total, residuals


