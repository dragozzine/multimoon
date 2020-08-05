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

"""
Inputs:
1) fit_array, the array of all the fitted parameters
Outputs:
1) log_likelihood, the log likelihood of the parameters with the priors
"""
def log_likelihood(params, obsdf, runprops, geo_obj_pos):
    # assuming Gaussian independent observations log-likelihood = -1/2 * chisquare
    
   # print(params, obsdf, geo_obj_pos)
    lh,residuals = mm_chisquare(params,obsdf, runprops, geo_obj_pos)
    lh = lh*-0.5
    #print('lh ',lh)

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
    objname = runprops.get("objectname")
    priorFilename = "../data/" +objname + "/" + objname + "_priors_df.csv"
    priors = pd.read_csv(priorFilename, sep='\t',index_col=0)
    priors = priors.transpose()
    
    name_dict = runprops.get("names_dict")
    
    
    params = mm_param.from_fit_array_to_param_df(float_params, float_names, fixed_df, total_df_names, fit_scale, name_dict)

    #print(params)
    
    lp = prior.mm_priors(priors,params,runprops)
    if runprops.get('verbose'):
        print('LogPriors: ',lp)

    if not np.isfinite(lp):
        return -np.inf
    #print(params)
    #print(obsdf)
    #print(geo_obj_pos)
    log_likeli, residuals = log_likelihood(params, obsdf, runprops, geo_obj_pos)
    llhood = lp + log_likeli

    if llhood > best_llhoods.get("best_llhood") and runprops.get("is_mcmc") and runprops.get("updatebestfitfile") :
        if runprops.get('verbose'):
            print("Previous best_llhoods, new llhood: ", best_llhoods.get('best_llhood'), llhood)
        best_llhoods['best_llhood'] = llhood
        best_llhoods['best_params'] = params.to_dict()
        the_file = runprops.get('runs_folder') + '/best_likelihoods.csv'

        reduced_chi_sq = llhood/(-0.5)/best_llhoods.get('deg_freedom')

        with open(the_file, 'a+', newline='') as write_obj:
            csv_writer = writer(write_obj, delimiter = ',')
            thelist = params.head(1).values.tolist()[0]
            thelist.insert(0, lp)
            thelist.insert(0, reduced_chi_sq)
            thelist.insert(0, llhood)
            for i in range(runprops.get('numobjects')):
                thelist.pop()
            for i in range(runprops.get("numobjects")-1):
                thelist.append(residuals[i])
                thelist.append(residuals[i+1])
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

    # Sorts them into ascending order
    import logging 
    try:
        time_arr_sec = time_arr*86400
        vec_df = generate_vector(paramdf, time_arr_sec, runprops)
    except:
        print('There was an error thrown within spinny:\n')
        logging.exception('')
        return np.inf
    names_dict = runprops.get("names_dict")
    names=[0 for i in range(numObj)]
    for i in range(0,numObj):
        names[i] = names_dict.get("name_"+str(i+1))
        
        
    # vec_df is a dataframe with len(time_arr) rows and 
    # columns are state parameters x nobjects
    # Example: vecdf["X_Pos_"+paramsdf["name_2"]] gets the x position of object 2
    # ecliptic (J2000) coordinates
    # km, kg, rad, s
    # primaricentric 

    name_1 = "X_Pos_"+names[0]
    
    if (vec_df[name_1][0] != 0.0):
        print("Not primaricentric like I thought!")
    
    Model_DeltaLong = np.zeros((numObj-1,len(time_arr)))
    Model_DeltaLat = np.zeros((numObj-1,len(time_arr)))
    #print(vec_df)

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

    for i in range(rows):
        for j in range(1,numObj):
#            #Check to make sure that these column names exist in the obsdf
#            if not names[j] in Model_DeltaLong.columns:
#                print(names[j], " is missing from the DeltaLong dataframe. Aborting run.")
#                print(Model_DeltaLong)
#            elif not "DeltaLong_"+names[j] in obsdf.columns:
#                print("DeltaLong_",names[j], " is missing from the DeltaLong dataframe. Aborting run.")
#                print(obsdf)
#            elif not "DeltaLong_"+names[j]+"_err" in obsdf.columns:
#                print("DeltaLong_",names[j], "_err is missing from the obsdf dataframe. Aborting run.")
#                print(obsdf)
#            elif not names[j] in Model_DeltaLat.columns: 
#                print(names[j], " is missing from the DeltaLat dataframe. Aborting run.")
#                print(Model_DeltaLat)
#            elif not "DeltaLat_"+names[j] in obsdf.columns:
#                print("DeltaLat_",names[j], " is missing from the obs dataframe. Aborting run.")
#                print(obsdf)
#            elif not "DeltaLat_"+names[j]+"_err" in obsdf.columns:
#                print("DeltaLat_",names[j], "_err is missing from the obs dataframe. Aborting run.")
#                print(obsdf)
#                sys.exit()
            
            #print('DLon', Model_DeltaLong[j-1][i],'\n obs_Dlon', obsdf["DeltaLong_"+names[j]][i],'\nobs_DLon_Err', obsdf["DeltaLong_"+names[j]+"_err"][i])
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


