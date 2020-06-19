from spinny_generate import *
from spinny_nosun import *
from keplerian import *
from spinny_plots import *
import numpy as np
import time
from time import ctime
import pandas as pd
import sys
sys.path.insert(1,"..")
import mm_runprops
runprops = mm_runprops.runprops

def generate_vector(paramsdf, t_arr):
    global G
    G = 6.674e-20 # Gravitational constant in km
    
    sys_df = paramsdf
    runprops = mm_runprops.runprops
    
    tol = runprops.get("spinny_tolerance")
    includesun = runprops.get("includesun")
        
    N = runprops.get("numobjects") # total number of objects in the system
    T = len(t_arr)                 # number of observation times

    j2s = []
    for col in sys_df.columns:
        if 'j2r2' in col:
            j2s.append(sys_df[col].iloc[0])
    
    j2_sum = sum(j2s)

    masses = []
    for col in sys_df.columns:
        if 'mass' in col:
            masses.append(sys_df[col].iloc[0])
    mass_sum = sum(masses)
    mass_primary = sys_df["mass_1"].iloc[0]
    
    if N == 2 and j2_sum == 0.00:  # checks if all objects are point masses, does keplerian integration instead
        kepler_system = kepler_2body(sys_df,t_arr)
        s_df = kepler_system[0]
        names = kepler_system[1]

    if mass_sum == mass_primary: #checks if only the primary has mass, all other bodies are massless
        kepler_system = kepler_Nbody(sys_df,t_arr)
        s_df = kepler_system[0]
        names = kepler_system[1]
        
    elif not "name_0" in sys_df.columns and not includesun: # runs a SPINNY integration without the sun if not included  
        system = build_spinny_multimoon(sys_df)
        spinny = evolve_spinny_ns(system[0],system[1],system[2],system[3],system[4],system[5],t_arr,tol)
        s_df = spinny[0]
        names = spinny[2]
        
    else:                         # runs SPINNY with the sun included
        system = build_spinny_multimoon(sys_df)
        spinny = evolve_spinny(system[0],system[1],system[2],system[3],system[4],system[5],t_arr)
        s_df = spinny[0]
        names = spinny[2]
        

    # creates a new dataframe using x,y,z position and velocity for each body
    data = {"Times":t_arr}
    
    for name in names:
        data.setdefault("X_Pos_"+name, s_df["X_Pos_"+name])
        data.setdefault("Y_Pos_"+name, s_df["Y_Pos_"+name])
        data.setdefault("Z_Pos_"+name, s_df["Z_Pos_"+name])
        data.setdefault("X_Vel_"+name, s_df["X_Vel_"+name])
        data.setdefault("Y_Vel_"+name, s_df["Y_Vel_"+name])
        data.setdefault("Z_Vel_"+name, s_df["Z_Vel_"+name])
        
    vec_df = pd.DataFrame(data)
  
    return(vec_df)
            

def build_spinny_multimoon(sys_df): 
    verbose = runprops.get("verbose")
    if verbose:
        print("Reading file to dataframe...")
        
    sys_df.reset_index()
    
    # in DESCENDING ORDER of size, an array of the masses 
    masses = sorted([sys_df[col].iloc[0] for col in sys_df.columns if "mass_" in col],reverse=True) 

    N = runprops.get("numobjects") #len(masses) # number of bodies in the system
    tol = runprops.get("spinny_tolerance")

    cols = []
    num = 0
    for col in sys_df.columns:
        if 'mass' in col:
            cols.append(num)
        num += 1

    
    #NOTE: as long as the runprops dict is sorted in order of descending mass, this should not be needed
    #cols = list(sys_df.columns[[int(np.where(sys_df==m)[1].flatten()) for m in masses]])
    #body_idx = [int(body[-1]) for body in cols]

    names_arr = np.empty(N,dtype="object")
    phys_arr = np.zeros((N,4))
    orb_arr = np.zeros((N,6))
    spin_arr = np.zeros((N,4))
    quat_arr = np.zeros((N,4))
    
    i = 0

    for n in range(1,N+1): # for each body in the system, added in order of descending mass:

        if "name_"+str(n) in sys_df.columns:
            names_arr[i] = sys_df["name_"+str(n)].iloc[0] # set name of the body
            
        else:
            names_arr[i] = "Body_"+str(n) 
        
        if "mass_"+str(n) in sys_df.columns:
            mass_n = sys_df["mass_"+str(n)].iloc[0] # set mass of the body
        else:
            mass_n = 0.0                    # default value--set to 0.0 if none given 
            
        if "ax_"+str(n) in sys_df.columns:
            ax_n = sys_df["ax_"+str(n)].iloc[0]
        else:
            ax_n = 1.0
            
        if "j2r2_"+str(n) in sys_df.columns:
            j2r2_n = sys_df["j2r2_"+str(n)].iloc[0]
        else:
            j2r2_n = 0.0    
        
        if "c22r2_"+str(n) in sys_df.columns:
            c22r2_n = sys_df["c22r2_"+str(n)].iloc[0]
        else:
            c22r2_n = 0.0    
                
        # set physical properties array
        phys_arr[i] = np.array([mass_n, ax_n, j2r2_n, c22r2_n])

        
        if "sma_"+str(n) in sys_df.columns:
            sma_n = sys_df["sma_"+str(n)].iloc[0]
        else:
            sma_n = 0.0
        
        if "ecc_"+str(n) in sys_df.columns:
            ecc_n = sys_df["ecc_"+str(n)].iloc[0]
        else:
            ecc_n = 0.0
            
        if "aop_"+str(n) in sys_df.columns:
             # convert all degree arguments to radians for SPINNY
            aop_n = (2*np.pi/180.)*sys_df["aop_"+str(n)].iloc[0]
        else:
            aop_n = 0.0
            
        if "inc_"+str(n) in sys_df.columns:
            inc_n = (np.pi/180.)*sys_df["inc_"+str(n)].iloc[0]
        else:
            inc_n = 0.0 
            
        if "lan_"+str(n) in sys_df.columns:
            lan_n = (np.pi/180.)*sys_df["lan_"+str(n)].iloc[0]
        else:
            lan_n = 0.0
        
        if "mea_"+str(n) in sys_df.columns:
            mea_n = (np.pi/180.)*sys_df["mea_"+str(n)].iloc[0]
        else:
            mea_n = 0.0
            
        # set orbital properties array
        orb_arr[i] = np.array([sma_n, ecc_n, aop_n, inc_n, lan_n, mea_n]) 
        
        i = i + 1
    # set orbit of primary equal to orbit of secondary (it will be scaled later)
    if "Sun" in names_arr:
        orb_arr[1] = orb_arr[2]
    else:
        orb_arr[0] = orb_arr[1]
    
    i = 0
    for n in range(1,N+1):
        # precession, obliquity angles are measured with respect to the ECLIPTIC, not the body's orbit
        # default values are set to be aligned with the orbit (LAN for prec, inc for obliq, AOP for longitude)
        #n= n-1
        if "splan_"+str(n) in sys_df.columns:
            sp_prc_n = (np.pi/180.)*sys_df["splan_"+str(n)].iloc[0]
        else:
            sp_prc_n = orb_arr[i,4]
            
        if "spinc_"+str(n) in sys_df.columns:
            sp_obl_n = (np.pi/180.)*sys_df["spinc_"+str(n)].iloc[0]
        else:
            sp_obl_n = orb_arr[i,3]
            
        if "spaop_"+str(n) in sys_df.columns:
            sp_lon_n = (np.pi/180.)*sys_df["spaop_"+str(n)].iloc[0]
        else:
            sp_lon_n = orb_arr[i,2]
            
        if "sprate_"+str(n) in sys_df.columns:
            sp_rate_n = sys_df["sprate_"+str(n)].iloc[0]
        else:
            # the default is equal to the orbital period, calculated from semi-major axis, sum of masses of body and primary 
            sp_rate_n = np.sqrt(orb_arr[i,0]**3.0/(G*(phys_arr[i,0]+phys_arr[1,0])) ) 
           
        # set spin properties array
        spin_arr[i] = np.array([sp_prc_n, sp_obl_n, sp_lon_n, sp_rate_n])

        # create orientation quaternions for each body from spin/orbit data
        
        # generate quaternion for obliquity angle of rotation axis
        # Euler angles: precession, obliquity, longitude, relative to the ecliptic
        
        r = R.from_euler('ZXZ', [spin_arr[i,0], spin_arr[i,1], spin_arr[i,2]])
        quat_i = r.as_quat()
        
        qi = quat_i[0]
        qj = quat_i[1]
        qk = quat_i[2]
        qr = quat_i[3]
        # set quaternion array
        # the order of the quaternion must be changed, since SPINNY needs scalar-first form
        quat_arr[i] = np.array([qr,qi,qj,qk]) #obliq_quat_n[0]

        i = i+1
        
    return(N, names_arr, phys_arr, orb_arr, spin_arr, quat_arr)


            
"""
spinny_output takes paramaters from the fitter and generates dataframe and plots of orbits, orbital parameters, and spin parameters over a specified period of time

The function returns a dataframe with all the bodies' orbital and spin parameters at each time step

## Future to do:
    - Flags passed through to specify which plots will be generated
    - Overplotting of multiple system models at once
    - Possibly add ability to create 3D VPython models?
    - Print out fractional change in energy/momentum if verbose flag?
""" 
def spinny_output(sys_df,t_start,t_end):
    
    runprops = mm_runprops.runprops
    
    t_arr = np.linspace(t_start,t_end,2000) # start and end times should be in SECONDS after a specified epoch
                                            # TODO: add code to adjust sampling rate for any given system   
    tol = runprops.get("spinny_tolerance")
    includesun = runprops.get("includesun")
        
    N = runprops.get("numobjects") # total number of objects in the system
    T = len(t_arr)                 # number of observation times

    j2s = []
    for col in sys_df.columns:
        if 'j2r2' in col:
            j2s.append(sys_df[col].iloc[0])
    
    j2_sum = sum(j2s)
    
    masses = []
    for col in sys_df.columns:
        if 'mass' in col:
            masses.append(sys_df[col].iloc[0])
    mass_sum = sum(masses)
    mass_primary = sys_df["mass_1"].iloc[0]
    name_primary = sys_df["name_1"].iloc[0]
    
    if N == 2 and j2_sum == 0.00:  # checks if all objects are point masses, does keplerian integration instead
        kepler_system = kepler_2body(sys_df,t_arr)
        s_df = kepler_system[0]
        names = kepler_system[1]

    if mass_sum == mass_primary: #checks if only the primary has mass, all other bodies are massless
        kepler_system = kepler_Nbody(sys_df,t_arr)
        s_df = kepler_system[0]
        names = kepler_system[1]
        
    elif not "name_0" in sys_df.columns and not includesun: # runs a SPINNY integration without the sun if not included  
        system = build_spinny_multimoon(sys_df)
        spinny = evolve_spinny_ns(system[0],system[1],system[2],system[3],system[4],system[5],t_arr,tol)
        s_df = spinny[0]
        names = spinny[2]
        
    else:                         # runs SPINNY with the sun included
        system = build_spinny_multimoon(sys_df)
        spinny = evolve_spinny(system[0],system[1],system[2],system[3],system[4],system[5],t_arr)
        s_df = spinny[0]
        names = spinny[2]        
    
    print("Generating .csv...")
    t_current = ctime().replace(" ","_")
    filename = name_primary+"_SPINNY_"+t_current+".csv"
    s_df.to_csv("../results/SPINNY-models/"+filename) ### I don't know if this path works--should save to results folder
    
    print("\n Generating figures...")
    spinny_plot(s_df, names)

    return(s_df)

def spinny_output_multiple(sys_df,t_start,t_end): #sys_df should have multiple rows with different parameters
    
    sys_df = pd.read_csv("../data/spinny_test_data/mm_multiple_test.csv")
    t_start = 0.
    t_end = 60.*24.*3600.0
    num_rows = len(sys_df.index) # total number of rows/models to generate
    t_arr = np.linspace(t_start,t_end,2000)
    tol = runprops.get("spinny_tolerance")
    N = runprops.get("numobjects")
    T = len(t_arr)     
     
    
    df_list = ["0"]*num_rows
    
    for i in range(0,num_rows):
        i_df = sys_df.iloc[[i]]
        
        system = build_spinny_multimoon(i_df)
        spinny = evolve_spinny_ns(system[0],system[1],system[2],system[3],system[4],system[5],t_arr,tol)
        s_df = spinny[0]
        names = spinny[2]
        
        df_list[i] = s_df
        
    spinny_plot_multiple(df_list,names,t_arr)
  
    
    
    
    
    

        