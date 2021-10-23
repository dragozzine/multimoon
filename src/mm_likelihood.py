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
import mpmath as mp
#from func_timeout import func_timeout, FunctionTimedOut

"""
Inputs:
1) fit_array, the array of all the fitted parameters
Outputs:
1) log_likelihood, the log likelihood of the parameters with the priors
"""
def log_likelihood(params, obsdf, runprops, geo_obj_pos):
    # assuming Gaussian independent observations log-likelihood = -1/2 * chisquare
    
    #print(params, obsdf, geo_obj_pos)
    lh,residuals = mm_chisquare(params,obsdf, runprops, geo_obj_pos)
    lh = lh*-0.5
    #print('lh ',lh)

    if runprops.get("robust_stats"):
        rows = obsdf.shape[0]
        numObj = runprops.get("numobjects")
        lh_robust = 0
        jitter = params["jitter"].iloc[0]
        p_outlier = params["pbad"].iloc[0]

        names_dict = runprops.get("names_dict")
        names=[0 for i in range(numObj)]
        for i in range(0,numObj):
            names[i] = names_dict.get("name_"+str(i+1))
        lh_robust_lon = 0
        lh_robust_lat = 0

        for j in range(1,numObj):
            for i in range(rows):
                combinedlon_err = np.sqrt(obsdf["DeltaLong_"+names[j]+"_err"][i]**2 + jitter**2)
                combinedlat_err = np.sqrt(obsdf["DeltaLat_"+names[j]+"_err"][i]**2 + jitter**2)
                omc_lon = (residuals[2*(j-1)][i] * obsdf["DeltaLong_"+names[j]+"_err"][i])**2
                omc_lat = (residuals[2*(j-1)+1][i] * obsdf["DeltaLat_"+names[j]+"_err"][i])**2
                #print(omc_lon,omc_lat)
                lh_robust_lon = mp.log( ((1-p_outlier)/(np.sqrt(2*np.pi*obsdf["DeltaLong_"+names[j]+"_err"][i]**2)))*mp.exp(-omc_lon/(2*obsdf["DeltaLong_"+names[j]+"_err"][i]**2)) + (p_outlier/np.sqrt(2*np.pi*combinedlon_err**2))*mp.exp(-omc_lon/(2*combinedlon_err**2))  )
                lh_robust_lat = mp.log( ((1-p_outlier)/(np.sqrt(2*np.pi*obsdf["DeltaLat_"+names[j]+"_err"][i]**2)))*mp.exp(-omc_lat/(2*obsdf["DeltaLat_"+names[j]+"_err"][i]**2)) + (p_outlier/np.sqrt(2*np.pi*combinedlat_err**2))*mp.exp(-omc_lat/(2*combinedlat_err**2))  )
                #print(names[j],lh_robust_lat,lh_robust_lon)

                if not (mp.isnan(lh_robust_lon) and mp.isnan(lh_robust_lat)):
                    #print(((1-p_outlier)/(np.sqrt(2*np.pi*obsdf["DeltaLong_"+names[j]+"_err"][i]**2)))*np.exp(-omc_lon/(2*obsdf["DeltaLong_"+names[j]+"_err"][i]**2)) + (p_outlier/np.sqrt(2*np.pi*combinedlon_err**2))*np.exp(-omc_lon/(2*combinedlon_err**2)))
                    #print(names[j],lh_robust_lat,lh_robust_lon)
                    lh_robust += lh_robust_lat + lh_robust_lon

        #print(lh_robust, lh, lh_robust-lh)
        return lh_robust, residuals

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
    
    #print('float_params read in from p0: \n',float_params)
    #print(float_params, float_names)
    objname = runprops.get("objectname")

    priorFilename = runprops.get('priors_filename')

    
    if 'runs' in os.getcwd() or 'results' in os.getcwd():
        os.chdir('../../../src')


    priors = pd.read_csv(priorFilename, sep='\t',index_col=0)
    priors = priors.transpose()
    
    name_dict = runprops.get("names_dict")
    if runprops.get('includesun') == True:
        name_dict['name_0'] = 'Sun'
    #print("floats:", float_params)
    #print("fixed:", fixed_df)
    #print("fit_scale:", fit_scale)
    params,fit_params = mm_param.from_fit_array_to_param_df(float_params, float_names, fixed_df, total_df_names, fit_scale, name_dict, runprops)
    #print('likelihood params 62', params)
    #if runprops.get('includesun') == 1:
    #    params.insert(0,'name_0',['Sun'])
        
    #print('Params: ',params)
    #print('Priors: ',fit_params)
    
    lp = prior.mm_priors(priors,params,runprops)
    #print('LP: ', lp)
    if runprops.get('verbose'):
        print('LogPriors: ',lp)
    if not np.isfinite(lp):
        return -np.inf

    log_likeli, residuals = log_likelihood(params, obsdf, runprops, geo_obj_pos)
    llhood = lp + log_likeli
    #print(llhood)
    the_file = runprops.get('results_folder') + '/best_likelihoods.csv'

    #You will notice this differes from the regular runs way to save data
    #Since we are using mpi, we need to continually retrieve the best_likelihoods csv
    #I did the math with some intense testing, and found this will only slow down
    #a 1000 step system by 1 minute, which typically takes 2 hours, so there is not much slow down.
    
    #print(llhood, best_llhoods.get('best_llhood'), runprops.get("is_mcmc"), runprops.get("updatebestfitfile"))
    if llhood > best_llhoods.get("best_llhood") and runprops.get("is_mcmc") and runprops.get("updatebestfitfile") :
        #print('is_mcmc')
        #if runprops.get('verbose'):
        #print("Previous best_llhoods, new llhood: ", best_llhoods.get('best_llhood'), llhood)
        best_llhoods['best_llhood'] = llhood
        #curr_best = best_llhoods["best_llhood"]
        #print(best_llhoods.get('best_llhood'))
        best_llhoods['best_params'] = params.to_dict()
        best_csv = pd.read_csv(the_file, index_col=None)        
        
        if len(best_csv.index) < 1:
            curr_best = -np.inf
        else:
            curr_best = best_csv.iloc[-1,0]
            #print('Curr_best:', curr_best)
            #print('Llhood:', llhood)
            
        #num_rows = len(best_csv.index)+1
        #print("Params: ", params)
        #print("Fit Params: ", fit_params)
        #print(llhood, curr_best)
        if llhood > curr_best:
        #if True:
            #print('adding')
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
                #thelist.insert(0, '')
                #print(thelist)
                for i in range(runprops.get('numobjects')):
                    thelist.pop()
                if runprops.get('includesun'):
                    #print('likelihood lin e130')
                    thelist.pop()
                #print(fit_params.head(1).values.tolist()[0])
                #print(fit_params)
                

                #for i in fit_params.head(1).values.tolist()[0]:
                #    thelist.append(i)


                for i in range(runprops.get("numobjects")-1):
                    thelist.append(residuals[2*(i-1)])
                    thelist.append(residuals[2*(i-1)+1])
                csv_writer.writerow(thelist)
                #print(thelist)
    #print(llhood)
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

    #print(paramdf)
    numObj = runprops.get("numobjects")
    verbose = runprops.get("verbose")
    pd.set_option('display.max_columns', None)
    names = []
    for i in range(1,numObj+1):
        names.append('name_'+str(i))
        if not 'name_'+str(i) in paramdf.columns:
            print('The parameter name_' + str(i)+ ' is not found in the parameter dataframe.')
            sys.exit()
        
    
    # use parameters dataframe with one set of parameters and observation times to call SPINNY to get the model dataframe
    # SPINNY returns (full?) state vector of each object (in primaricentric ecliptic J2000 coordinates) at each time
    
    # parameters dataframe "paramdf" is 1 row and many columns
    # N objects, each of which is a point-mass, oblate, or triaxial
       # point-masses: mass, 6 orbital elements
       # oblate: + J2R2, spinc, splan
       # triaxial: + C22R2, spaop, sprate
       # dynamicstoincludeflags='2010'
    # Include Sun or not
    # it has already been "checked" to make sure that all of the columns that are needed are there
    # mass_0, sma_0, ecc_0, ..., mea_0 are the orbital elements of the Sun 
       # at the epoch in primaricentric coordinates
    # when columns are not needed they are typically absent

    # MultiMoon units are km, kg, deg, seconds

    obsdf = obsdf.sort_values(by=['time'])
    #time_arr = np.sort(obsdf['time'].values.flatten())
    time_arr = obsdf['time'].values.flatten()# gets an array of observation times from the obs dataframe

    # Setting times relative to the epoch
    epoch = runprops.get("epoch_SJD")
    time_arr = time_arr - epoch

    # Sorts them into ascending order
#    import logging
    #print(paramdf)

    begin = time.time()
    #print(paramdf)
    try:
        time_arr_sec = time_arr*86400
        #vec_df = func_timeout(5,generate_vector,args=(paramdf, time_arr_sec, runprops))
        #print(paramdf)
        vec_df = generate_vector(paramdf, time_arr_sec, runprops)
    #except FunctionTimedOut:
    #    print('Spinny took longer than 5 seconds to run 1 walker-step:\n')
    #    return np.inf
#    except Exception as e:
#        logging.exception('')
#        return np.inf
    except Exception as e:
        print('There was an error thrown within spinny:\n', e)
        return -np.inf
    names_dict = runprops.get("names_dict")
    names=[0 for i in range(numObj)]
    for i in range(0,numObj):
        names[i] = names_dict.get("name_"+str(i+1))
    end = time.time()
    #print(end-begin,' seconds')
    #if end-begin > 2:
    #    print('A step took ', end-begin,' seconds')
    #    print(paramdf)

    # vec_df is a dataframe with len(time_arr) rows and
    # columns are state parameters x nobjects
    # Example: vecdf["X_Pos_"+paramsdf["name_2"]] gets the x position of object 2
    # ecliptic (J2000) coordinates
    # km, kg, rad, s
    # primaricentric 

    name_1 = "X_Pos_"+names[0]
    #print(vec_df)
    if (vec_df[name_1][0] != 0.0):
        print("Not primaricentric like I thought!")
        #print("vec_df[name_1] = ", vec_df)
    
    Model_DeltaLong = np.zeros((numObj-1,len(time_arr)))
    Model_DeltaLat = np.zeros((numObj-1,len(time_arr)))
    if runprops.get('includesun') == 1:
        #print(vec_df)
        vec_df = vec_df.drop(['X_Pos_Sun', 'Y_Pos_Sun', 'Z_Pos_Sun', 'X_Vel_Sun', 'Y_Vel_Sun', 'Z_Vel_Sun'], axis=1)

    positionData = np.zeros((numObj*3,len(time_arr)))
        
    for i in range(0,numObj):
        positionData[3*i] = vec_df["X_Pos_"+names[i]]
        positionData[3*i+1] = vec_df["Y_Pos_"+names[i]]
        positionData[3*i+2] = vec_df["Z_Pos_"+names[i]]

        # tind = index/row number of vec_df corresponding to this time
        
        # for object from 2:N
             # gather relative positions
             # thisdx = vec_df[tind,"X_Pos_"+paramdf["name_"+str(object)] [ - X_Pos_1=0]
             # same for dy and dz

        # Implicitly assume that observer is close to geocenter (within ~0.01 AU)
        
        # obs_to_prim_pos = vectors of observer to primary
        # prim_to_sat__pos = vectors of primary to satellite

    obs_to_prim_pos = [positionData[0]+geo_obj_pos['x'].tolist(),positionData[1]+geo_obj_pos['y'].tolist(),positionData[2]+geo_obj_pos['z'].tolist()]
    for i in range(1,numObj):
        prim_to_sat_pos = [positionData[i*3],positionData[i*3+1],positionData[i*3+2]]
        Model_DeltaLong[i-1], Model_DeltaLat[i-1] = mm_relast.convert_ecl_rel_pos_to_geo_rel_ast(obs_to_prim_pos, prim_to_sat_pos)
        

        if verbose: 
            print('i', i)
            print('obs_to_prim', obs_to_prim_pos,'\nprim_to', prim_to_sat_pos)
            print('DeltaLong', Model_DeltaLong, '\nDeltaLat', Model_DeltaLat)
        # mm_relast
        
        # obs_to_prim_pos = vector position of the observer relative to the primary (in J2000 ecliptic frame)
          # input: numpy array of 3 values: x, y, z in units of km
          # comes from geocentric_object_position dataframe (defined in mm_run)
          # DS TODO: change convert_ecl_rel_pos... to input in km and change geo_object to be in km too
                  
        # prim_to_sat_pos = position of the primary relative to the satellite (in J2000 ecliptic)
          # numpy array of 3 values: x, y, z in units of km
          # comes from thisdx, thisdy, thisdz from vec_df above

        # Output: Delta Long, Delta Lat in arcseconds
        # Here is where you would correct for photocenter shifts
        # model delta long and delta lat
        # Put into model dataframe? 

    # Putting in COM-COL offset
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

        
    
    
    if runprops.get("com_offset"):
        Model_DeltaLong = Model_DeltaLong + paramdf["long_offset"].iloc[0]
        Model_DeltaLat = Model_DeltaLat + paramdf["lat_offset"].iloc[0]

    # Outputting the Model_DeltaLong and Lat if gensynth flag is included in function call
    # This was done by BP. PLEASE let me know if this breaks things
    if gensynth:
        if verbose:
            print("Returning the Model_DeltaLong and Lat dataframes for use in synthetic astrometry.")
        return Model_DeltaLong, Model_DeltaLat, obsdf

    # Now we have model delta Long and delta Lat for each object and each time 
    rows = obsdf.shape[0]

    residuals = np.zeros(((numObj-1)*2, rows))
    get_residuals = runprops.get("get_resid")
    delta_offset = 0
    #if runprops.get('search_inner'):
        
        
        
    for i in range(rows):
        for j in range(1,numObj):

            residuals[2*(j-1)][i] = ((Model_DeltaLong[j-1][i]-obsdf["DeltaLong_"+names[j]][i])/obsdf["DeltaLong_"+names[j]+"_err"][i])
            residuals[2*(j-1)+1][i] = ((Model_DeltaLat[j-1][i]-obsdf["DeltaLat_"+names[j]][i])/obsdf["DeltaLat_"+names[j]+"_err"][i])
           
                
            if verbose:
                print("i,j,model,obs,err")
                print(i, j, Model_DeltaLong[j-1][i], obsdf["DeltaLong_"+names[j]][i], obsdf["DeltaLong_"+names[j]+"_err"][i])

                                      
    # Loop through obsdf and for each defined value of delta Long/Lat 
    # calculate chisquare = sum [ (model-obs)/err ] ^2 
    # associate the correct satellite in the observations with the objects from the model
    # using the name of the object "Model_DeltaLong_Hiiaka" with "DeltaLong_Hiiaka" 
    # and "DeltaLong_Hiiaka_err"
    # AND throw an error if the names don't all line up right
    
    
    chisquares = residuals**2
    #print('residuals ',residuals)
    
    chisq_tot = np.zeros(2*numObj)
    for i in range(0,2*numObj-2):
        chisq_tot[i]=np.nansum(chisquares[i])
        
    chisquare_total = np.nansum(chisq_tot)

    if verbose:
        print("chisq_tot, chisquare_total, residuals")
        print(chisq_tot, chisquare_total, residuals)

    # return chisquare

    
    return chisquare_total, residuals


