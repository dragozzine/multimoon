import pandas as pd
import numpy as np
import scipy as sp
import time

'''
    NAME:
         mm_priors
         
    PURPOSE:
         Given a dataframe of priors, and a dataframe of observed parameters for the same data, 
         this function calculates the likelihood of the distribution.
         
    CALLING SEQUENCE:
         totalLogProb = mm_priors(priors, params)
   
    INPUTS
         priors - a dataframe of 9 rows which holds the prior for the data points given, and the distribution shape wanted.
         params - A single rowed dataframe of the actual observed parameters for the object. 
          
    OUTPUTS:
          totalLogProb - The total Log of the probability of all the priors against the parameters
'''
def mm_priors(priors, params, runprops):
    columnList = list(priors)
    totalLogProb = 0

    # Remove any rows of priors that arent in params
    columnList = list(params)
    if "name_0" in columnList:
        columnList.remove("name_0")
    if "name_1" in columnList:
        columnList.remove("name_1")
    if "name_2" in columnList:
        columnList.remove("name_2")
    if "name_3" in columnList:
        columnList.remove("name_3")
    if "name_4" in columnList:
        columnList.remove("name_4")
    if "name_5" in columnList:
        columnList.remove("name_5")
    
    #This loop is to make sure all of the column values are floats, because pandas sometimes turns the values to strings when read from file
    for i in columnList:
        priors[i].astype(float)

    count = 0
    allProbs = []
    numNaNs = 0

    #This loop runs through every column in the priors dataframe, and evaluates the probability density
    #function of the specified type.
    for i in columnList:
        count += 1
        theInt = int(priors[i][0])

        for j in range(runprops.get('numobjects')):
            #Uniform Distribution Shape
            if theInt == 0:
                    
                if '_'+str(j+1) in i and params[i][0] <= priors[i][2] and params[i][0] >= priors[i][1]:
                    allProbs.append(1)
                elif '_'+str(j+1) in i and np.isnan(params[i][0]):
                    numNaNs += 1
                elif '_'+str(j+1) in i:
                    allProbs.append(0)
                    #print(i, " is log 0", " because i is ", params[i][0])
            
            #Log-Uniform Distribution Shape
            elif theInt == 1:
                if '_'+str(j+1) in i and params[i][0] < priors[i][4] and params[i][0] > priors[i][3]:
                    allProbs.append(1)
                elif '_'+str(j+1) in i and np.isnan(params[i][0]):
                    numNaNs += 1
                elif '_'+str(j+1) in i:
                    allProbs.append(0)
                
            # Normal Distribution Shape
            elif theInt == 2:
                if '_'+str(j+1) in i and not np.isnan(params[i][0]):
                    allProbs.append(np.exp(-1/2*((params[i][0]-priors[i][6])/priors[i][5])**2))
                
            #Log Normal Distribution Shape
            elif theInt == 3:
                if '_'+str(j+1) in i and not np.isnan(params[i][0]):
                    allProbs.append(np.exp(-1/2*(((np.log(params[i][0])-priors[i][8])**2)/(priors[i][7])**2))/params[i][0])
            else:
                a = 1 #print('N/A input for column: ', i, ' in priors dataframe.') 
            
            #Make sure the values in the params df are real.
            
            if 'mass' in i:
                if i in params and params[i][0] < 0:
                    #print(i, " is outside of the realistic value with a value of ", params[i][0])
                    return -np.inf
            elif 'ecc' in i:
                if i in params and params[i][0] < 0:
                    #print(i, " is outside of the realistic value with a value of ", params[i][0])
                    return -np.inf
                elif i in params and params[i][0] > 1:
                    #print(i, " is outside of the realistic value with a value of ", params[i][0])
                    return -np.inf
            elif 'sma' in i:
                if i in params and params[i][0] < 0:
                    #print(i, " is outside of the realistic value with a value of ", params[i][0])
                    return -np.inf
            elif 'j2r2' in i:
                if i in params and params[i][0] < 0:
                    #print(i, " is outside of the realistic value with a value of ", params[i][0])
                    return -np.inf      
            elif 'c22r2' in i:
                if i in params and params[i][0] < 0:
                    #print(i, " is outside of the realistic value with a value of ", params[i][0])
                    return -np.inf      
            #Here, add the Prior Probability Density function for this element to the total

    # Checking robust statistics priors
    if runprops.get("robust_stats"):
        logjitter = params["logjitter"].iloc[0]
        p_outlier = params["pbad"].iloc[0]
        
        # log(jitter) must be between -3.0-0.0 (0.001-1.0 arcseconds)
        if logjitter < -3.0:
            return -np.inf
        if logjitter > 0.0:
            return -np.inf
        # p_outlier must be between 0-1
        if p_outlier < 0.0:
            return -np.inf
        if p_outlier > 1.0:
            return -np.inf
        
    # Making sure that c22 < 0.5*j2
    dynamicstoincludeflags = runprops.get("dynamicstoincludeflags")
    for i in range(runprops.get("numobjects")):
        if dynamicstoincludeflags[i] == "2":
            if (params["j2r2_" + str(i+1)].values[0]*0.5 < params["c22r2_" + str(i+1)].values[0]):
                return -np.inf

    # Checking object specific properties...
    min_periapse = runprops.get("min_periapse")

    # Check that objects are ordered correctly.
    for i in range(2,runprops.get('numobjects')):
        if params['sma_'+str(i)].values > params['sma_'+str(i+1)].values:
            print('Objects should be input from closest to furthest object in orbit. We detect that your satellites are not input in this order right now in your initial guess folder. Please change this before running again.')
            import sys
            sys.exit()
    
    # Ensure sat-spin inc is <90 (forces prograde orbits). This can be removed to loosen this restriction.
    # TODO: Make this a setting in runprops?
    for i in range(1,runprops.get("numobjects")):
        if dynamicstoincludeflags[0] == "1" or dynamicstoincludeflags[0] == "2":
            spinc_1_n = params["spinc_1"]/180*np.pi
            splan_1_n = params["splan_1"]/180*np.pi
            inc_i_n = params["inc_"+str(i+1)]/180*np.pi
            lan_i_n = params["lan_"+str(i+1)]/180*np.pi
            
            mutualinc = np.arccos( np.cos(spinc_1_n)*np.cos(inc_i_n) + np.sin(spinc_1_n)*np.sin(inc_i_n)*np.cos(splan_1_n - lan_i_n) )
            mutualinc = np.rad2deg(mutualinc).values
            if mutualinc > 90:
                print('mutualinc > 90')
                return -np.inf
        
        # Making sure min periapse is obeyed        
        if i == 1 and (params["sma_" + str(i+1)].values[0]*(1-params["ecc_" + str(i+1)].values[0]) < min_periapse):
            return -np.inf

        # Enforce max apoapse for when there are two moons. Ensures collisions are not possible.
        elif i != 1 and (params["sma_" + str(i+1)].values[0]*(1-params["ecc_" + str(i+1)].values[0])-params["sma_" + str(i)].values[0]*(1+params["ecc_" + str(i)].values[0]) < min_periapse):
            return -np.inf

    # Makes sure mass1 is greater than mass2
    # TODO: Is this necessary???
    for i in range(1,runprops.get("numobjects")):
        mass1 = params["mass_" + str(1)].values[0]
        mass2 = params["mass_" + str(i+1)].values[0]
        if mass1 < mass2:
            if runprops.get('is_mcmc') == False:
                print('It seems one of your initial guesses prodcues a mass1 < mass2. MultiMoon currently runs best when the most massive object is the primary. If your run will not start due to a maximum reset, we recommend modifying your initial guesses to have the most massive object exist as the primary.')
            return -np.inf
    
    if runprops.get('verbose'):
        print('AllProbs:' ,allProbs)
    
    for x in allProbs:
        totalLogProb = totalLogProb + np.log(x)
  
    return totalLogProb


