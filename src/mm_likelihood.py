import numpy as np
import mm_priors as prior
import pandas as pd
import mm_param

"""
Inputs:
1) fit_array, the array of all the fitted parameters

Outputs:
1) log_likelihood, the log likelihood of the parameters with the priors

"""
def log_likelihood(params, obsdf, runprops, fitarray_to_params_dict):
    # assuming Guassian independent observations log-likehood = -1/2 * chisqure

    # DS TODO: convert parameters array to parameters dataframe
    
    #We currently do not have a defined constraints dictionary anywhere, so currently we'll define it as empty
    constraints = {}
    
    paramdf = mm_param.from_fit_array_to_param_df(params, constraints, fitarray_to_params_dict)

    lh = mm_chisquare(paramdf,obsdf)*-0.5
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
def log_probability(params, runprops, fitarray_to_params_dict, obsdf):
    
    objname = runprops.get("objname")
    priors = pd.read_csv("../data/" +objname + "/" + objname + "_priors_df.csv", sep='\t',index_col=0)
    priors = priors.transpose()
    
    lp = prior.mm_priors(priors,params)
    
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(params, obsdf, runprops, fitarray_to_params_dict)


"""
Inputs:
1)The Parameters dataframe
2) The Observation Dataframe

Outputs:
1) The chi-squared number of the likelihood
"""
# calculates the chi-square for parameters given observations
def mm_chisquare(paramdf, obsdf):


    # DS TODO: make sure "name_i" is in paramdf    
    
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


    time_arr = obsdf['time'].values.flatten() # gets an array of observation times from the obs dataframe
    # SP TODO: should times be sorted?

    # SP TODO: make sure SPINNY can handle i/o in MultiMoon units
    # SP TODO: return full state vector
    vec_df = spinny_vector(paramdf, time_arr) # returns a dataframe of state vectors for each body in the system
    # vec_df is a dataframe with len(time_arr) row and 
    # columns are state parameters x nobjects
    # Example: vecdf["X_Pos_"+paramsdf["name_2"]] gets the x position of object 2
    # ecliptic (J2000) coordinates
    # km, kg, rad, s
    # primaricentric 

    if (vec_df[0,"X_Pos_"+paramdf["name_1"]] != 0.0):
        print("Not primaricentric like I thought!")

    
    for t in time_arr:

        # DS TODO: get relative positions out of vec_df

        # tind = index/row number of vec_df corresponding to this time

        # for object from 2:N
             # gather relative positions
             # thisdx = vec_df[tind,"X_Pos_"+paramdf["name_"+str(object)] [ - X_Pos_1=0]
             # same for dy and dz

        # which numbers do I want from the obsdf?



        # Implicitly assume that observer is close to geocenter (within ~0.01 AU)


        convert_ecl_rel_pos_to_geo_rel_ast(ecl_rel_pos, obj_rel_pos, rel_moon):
        # DS TODO: study and update the code and comments
        # mm_relast
        # ecl_rel_pos = position of the primary relative to the observer (in ecliptic frame)
          # input: numpy array of 3 values: x, y, z in units of km
          # comes from geocentric_object_position dataframe (defined in mm_run)
          # DS TODO: change convert_ecl_rel_pos... to input in km and change geo_object to be in km too
          
        # obj_rel_pos = relative positions between object and primary (in ecliptic)
          # numpy array of 3 values: x, y, z in units of km
          # comes from thisdx, thisdy, thisdz from vec_df above

        # rel_moon
          # 

        # Output: Delta Long, Delta Lat in arcseconds
       

        # Here is where you would correct for photocenter shifts


        # model delta long and delta lat
        # Put into model dataframe? 


    # DS TODO:
    # Now we have model delta Long and delta Lat for each object and each time 
   
    # Loop through obsdf and for each defined value of delta Long/Lat 
    # calculate chisquare = sum [ (model-obs)/err ] ^2 
    # associate the correct satellite in the observations with the objects from the model
    # using the name of the object "Model_DeltaLong_Hiiaka" with "DeltaLong_Hiiaka" 
    # and "DeltaLong_Hiiaka_err"
    # AND throw an error if the names don't all line up right


    # return chisquare
    
    return 1
