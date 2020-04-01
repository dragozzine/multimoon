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
def mm_priors(priors, params):
    columnList = list(priors)
    totalLogProb = 0
    
    probDist = pd.DataFrame(columns = ['mass_1','mass_2','sma_2','ecc_2','aop_2','inc_2','lan_2','mea_2','mass_3','sma_3','ecc_3','aop_3','inc_3','lan_3','mea_3','j2r2_1','c22r2_1','spaop_1','spinc_1','splan_1','sprate_1'],index=['PDF'])
   
    count = 0
    allProbs = []
    numNaNs = 0
    #This loop runs through every column in the priors dataframe, and evaluates the probability density
    #function of the specified type.
    for i in columnList:
        count += 1
        
        #Uniform Distribution Shape
        if priors[i][0] == 0:
            if params[i][0] < priors[i][2] and params[i][0] > priors[i][1]:
                allProbs.append(1)
            elif np.isnan(x[count]):
                numNaNs += 1
            else:
                allProbs.append(0)
        
        #Log-Uniform Distribution Shape
        elif priors[i][0] == 0:
            if params[i][0] < priors[i][4] and params[i][0] > priors[i][3]:
                allProbs.append(1)
            elif np.isnan(params[i][0]):
                numNaNs += 1
            else:
                allProbs.append(0)
            
        # Normal Distribution Shape
        elif priors[i][0] == 2:
            if not np.isnan(params[i][0]):
                allProbs.append(np.exp(-1/2*((params[i][0]-priors[i][6])/priors[i][5])**2))
            
        #Log Normal Distribution Shape
        elif priors[i][0] == 3:
            if not np.isnan(params[i][0]):
                allProbs.append(np.exp(-1/2*(((np.log(params[i][0])-priors[i][8])**2)/(priors[i][7])**2))/params[i][0])
        else:
            print('Invalid input for: ', i) 
            
        #Here, add the Prior Probability Density function for this element to the total
    for x in allProbs:
        totalLogProb = totalLogProb + np.log(x)
        print(x)
        
    return totalLogProb