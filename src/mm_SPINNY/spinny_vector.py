from spinny_generate import *
from spinny_nosun import *
from keplerian import kepler_integrate, kepler_plot
import numpy as np
import time
from time import ctime
import pandas as pd
import sys


def generate_vector(paramsdf, t_arr):
   
    # need to rearrange the paramsdf into a shape that is easier for SPINNY to read?
    # TODO: Update this for when I know what the params dataframe looks like
    sys_df = paramsdf
    
    N = len(sys_df.columns) # total number of objects in the system
    T = len(t_arr)          # number of observation times
    
    j2_sum = sum(sys_df.loc["j2r2",:].values.flatten())
    names = list(sys_df.columns)
    
    if N == 2 and j2_sum == 0.00:  # checks if all objects are point masses, does keplerian integration instead
        kepler_system = kepler_integrate(sys_df,t_arr)
        kepler_df = kepler_system[0]
        names = kepler_system[1]

    elif not "Sun" in names:       # runs a SPINNY integration without the sun if not included  
        system = build_spinny_ns(sys_df)
        spinny = evolve_spinny_ns(system[0],system[1],system[2],system[3],system[4],system[5],t_arr)
        s_df = spinny[0]
        names = spinny[2]
        
    else:                         # runs SPINNY with the sun included
        system = build_spinny(sys_df)
        spinny = evolve_spinny(system[0],system[1],system[2],system[3],system[4],system[5],t_arr)
        s_df = spinny[0]
        names = spinny[2]
        
    for n in range(0,N):
 
        # creates a new dataframe using just x,y,z position for each body    
        vec_df = s_df[["X_Pos_"+name,"Y_Pos_"+name,"Z_Pos_"+name for name in names]].copy() ## Does this work??

    return(vec_df)
            
            
            
            
            
            
            
            
        