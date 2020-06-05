import numpy as np
import mm_priors as prior
import pandas as pd
import mm_param
import sys
sys.path.insert(1, 'mm_SPINNY/')
from mm_SPINNY.spinny_vector import generate_vector
import random
import mm_relast

"""
Inputs:
1) fit_array, the array of all the fitted parameters

Outputs:
1) log_likelihood, the log likelihood of the parameters with the priors

"""
def log_likelihood(params, obsdf, runprops):
    # assuming Gaussian independent observations log-likelihood = -1/2 * chisquare

    if runprops.get("get_resid"):
        lh,residuals = mm_chisquare(params,obsdf, runprops)*-0.5
    else:
        lh = mm_chisquare(params,obsdf, runprops)*-0.5
    
    return lh


"""
Inputs:
1) params
2) runprops
3) fitarray_dict
4) 

Outputs:
1) log_probability, the log_likelihood plus the priors, which is the total probability

"""
def log_probability(float_params, float_names, fixed_df, total_df_names, fit_scale, runprops, obsdf):
    
    objname = runprops.get("objectname")
    priorFilename = "../data/" +objname + "/" + objname + "_priors_df.csv"
    priors = pd.read_csv(priorFilename, sep='\t',index_col=0)
    priors = priors.transpose()
    
    name_dict = runprops.get("names_dict")
    params = mm_param.from_fit_array_to_param_df(float_params, float_names, fixed_df, total_df_names, fit_scale, name_dict)
    
    lp = prior.mm_priors(priors,params)
    
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, obsdf, runprops)


"""
Inputs:
1)The Parameters dataframe
2) The Observation Dataframe

Outputs:
1) The chi-squared number of the likelihood
"""
# calculates the chi-square for parameters given observations
def mm_chisquare(paramdf, obsdf, runprops):

    numObj = runprops.get("numobjects")
    verbose = runprops.get("verbose")
    print(verbose, verbose == True)
    if verbose: 
        print("verbose test works")
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


    time_arr = np.sort(obsdf['time'].values.flatten()) # gets an array of observation times from the obs dataframe
                                                       # Sorts them into ascending order

    #vec_df = generate_vector(paramdf, time_arr)
    names=[0 for i in range(numObj)]
    for i in range(0,numObj):
        names[i] = paramdf["name_"+str(i+1)][0]

    vec_df = pd.DataFrame(index = range(len(time_arr)))
    for i in names:
        x = "X_Pos_"+i
        y = "Y_Pos_"+i
        z = "Z_Pos_"+i
        vec_df[x] = random.sample(range(50000), len(time_arr))
        vec_df[y] = random.sample(range(50000), len(time_arr))
        vec_df[z] = random.sample(range(50000), len(time_arr))
    
    vec_df['time'] = time_arr
    # vec_df is a dataframe with len(time_arr) rows and 
    # columns are state parameters x nobjects
    # Example: vecdf["X_Pos_"+paramsdf["name_2"]] gets the x position of object 2
    # ecliptic (J2000) coordinates
    # km, kg, rad, s
    # primaricentric 
    
<<<<<<< HEAD
    name_1 = "X_Pos_"+paramdf["name_1"][0]
    if (vec_df[name_1][0] != 0.0):
        print("Not primaricentric like I thought!")
=======
    # SP: NOTE I commented this out because it was causing a TypeError -- fix later
    #if (vec_df[0,"X_Pos_"+paramdf["name_1"]] != 0.0):
    #    print("Not primaricentric like I thought!")
>>>>>>> 40655afffc19d0ae7d8507e2f3778d5ac4ac8ed3

    Model_DeltaLong = pd.DataFrame(index = range(len(time_arr)),columns = names)
    Model_DeltaLat = pd.DataFrame(index = range(len(time_arr)), columns = names)
    Model_DeltaLong['time'] = time_arr
    Model_DeltaLat['time'] = time_arr
    
    for t in range(len(time_arr)):
        # DS TODO: get relative positions out of vec_df
<<<<<<< HEAD
        positionData = pd.DataFrame()
      
        for i in range(0,numObj):
            positionData[names[i]] = vec_df["X_Pos_"+names[i]]
            positionData[names[i]] = vec_df["Y_Pos_"+names[i]]
            positionData[names[i]] = vec_df["Z_Pos_"+names[i]]
                                                               
=======
        positionData = np.empty(numObj,3)
        names = np.empty(numObj)
        for i in range(1,numObj+1):
            names.append(paramdf["name_"+str(i)])
            positionData[i][0] = vec_df["X_Pos_"+names[i]]
            positionData[i][1] = vec_df["Y_Pos_"+names[i]]
            positionData[i][2] = vec_df["Z_Pos_"+names[i]]
                                      
                                      
>>>>>>> 40655afffc19d0ae7d8507e2f3778d5ac4ac8ed3
        # tind = index/row number of vec_df corresponding to this time
        
        # for object from 2:N
             # gather relative positions
             # thisdx = vec_df[tind,"X_Pos_"+paramdf["name_"+str(object)] [ - X_Pos_1=0]
             # same for dy and dz

        # which numbers do I want from the obsdf?

        # Implicitly assume that observer is close to geocenter (within ~0.01 AU)
        
        # obs_to_prim_pos = vectors of observer to primary
        # prim_to_sat__pos = vectors of primary to satellite
        
        obs_to_prim_pos = positionData[names[0]]
        for i in range(1,numObj):
            Model_DeltaLong[names[i]][t], Model_DeltaLat[names[i]][t] = mm_relast.convert_ecl_rel_pos_to_geo_rel_ast(obs_to_prim_pos, positionData[names[i]])
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


    # Now we have model delta Long and delta Lat for each object and each time 
    rows = len(obsdf.iloc[0])
    chisquare = pd.DataFrame(columns = names, index = range(len(time_arr)))
    residuals = pd.DataFrame(columns = names, index = range(len(time_arr)))
    get_residuals = runprops.get("get_resid")
    print('rows: ',rows)
    for i in range(rows):
        for j in range(1,numObj):
            #Check to make sure that these column names exist in the obsdf
            if not names[j] in Model_DeltaLong.columns:
                print(names[j], " is missing from the DeltaLong dataframe. Aborting run.")
                print(Model_DeltaLong)
            elif not "DeltaLong_"+names[j] in obsdf.columns:
                print("DeltaLong_",names[j], " is missing from the DeltaLong dataframe. Aborting run.")
                print(obsdf)
            elif not "DeltaLong_"+names[j]+"_err" in obsdf.columns:
                print("DeltaLong_",names[j], "_err is missing from the obsdf dataframe. Aborting run.")
                print(obsdf)
            elif not names[j] in Model_DeltaLat.columns: 
                print(names[j], " is missing from the DeltaLat dataframe. Aborting run.")
                print(Model_DeltaLat)
            elif not "DeltaLat_"+names[j] in obsdf.columns:
                print("DeltaLat_",names[j], " is missing from the obs dataframe. Aborting run.")
                print(obsdf)
            elif not "DeltaLat_"+names[j]+"_err" in obsdf.columns:
                print("DeltaLat_",names[j], "_err is missing from the obs dataframe. Aborting run.")
                print(obsdf)

                sys.exit()
            

            residuals[names[j]][i] = ((Model_DeltaLong[names[j]][i]-obsdf["DeltaLong_"+names[j]][i])/obsdf["DeltaLong_"+names[j]+"_err"][i])**2
            residuals[names[j]][i] = ((Model_DeltaLat[names[j]][i]-obsdf["DeltaLat_"+names[j]][i])/obsdf["DeltaLat_"+names[j]+"_err"][i])**2
            chisquare[names[j]][i] = ((Model_DeltaLong[names[j]][i]-obsdf["DeltaLong_"+names[j]][i])/obsdf["DeltaLong_"+names[j]+"_err"][i])**2
            chisquare[names[j]][i] = ((Model_DeltaLat[names[j]][i]-obsdf["DeltaLat_"+names[j]][i])/obsdf["DeltaLat_"+names[j]+"_err"][i])**2
        
                                      
    # Loop through obsdf and for each defined value of delta Long/Lat 
    # calculate chisquare = sum [ (model-obs)/err ] ^2 
    # associate the correct satellite in the observations with the objects from the model
    # using the name of the object "Model_DeltaLong_Hiiaka" with "DeltaLong_Hiiaka" 
    # and "DeltaLong_Hiiaka_err"
    # AND throw an error if the names don't all line up right
    
    chisq_tot = pd.DataFrame(columns = names, index = [0])
    for i in range(1,numObj):
        chisq_tot[names[i]]=chisquare[names[i]].sum()
        
    chisquare_total = chisq_tot.sum(axis = 1, skipna = True)[0]

    # return chisquare
    if get_residuals:
        return chisquare_total, residuals
    else:
        return chisquare_total
