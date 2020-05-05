# testsupport.py
#
# Darin Ragozzine
# March 27, 2020
# 
# Supports testing of MultiMoon by including functions
# that return valid "shapes" and types of various
# planned pieces of MultiMoon. 

import pandas
import numpy as np

# test_runprops
# 

def test_runprops():
    """ Returns a runprops dictionary. """
    
    runprops = {
        'date': 'today',
        'time': 'now',
        'user': 'autogenerate',
        'MultiMoon commit hash': 'autogenerate',
        'mm_initguess_which' : 'get_init_guess_from_dist', # which function from mm_initguess to use
        'numobjects' : 3,
        'RunGoal' : 'user puts their goal for this run here',
        'objectname' : 'Haumea',
        'nwalkers' : 30,
        'nsteps' : 1001,
        'nthinning' : 100
    }
    
#Filename for starting guess - default = start_guess.csv
#fixfloat_df - default="all float" dataframe in parameters format
#lpriorfunc - default="log_prior"
#llhoodfunc - default="log_likelihood"
#obsdata - name of file with text/csv of observations dataframe
#posteriorfile - name of file for posterior, default="date_time_objname.csv"
#verbose - default=1, how verbose to be in terms of output
#(0 = quiet, no text; 1 = typical, only major outputs, 5 = everything) 
#plotflags: True if this plot is to be made, False if not
##plot_model_in_data - show model in data space
#plot_model_over_time - show model continuous in time over data space
#other plots too

    return runprops


def test_mmparamdf():
    # Edited Seth Pincock 
    # March 30, 2020
    # Modified df to include parameters from Haumea system 
    
    """ Returns a test parameters dataframe with fake values."""

    objects=["0","1","2","3"] # 3 objects in KBO system and Sun (0)
    objparams=["name","mass","sma","ecc","aop","inc","lan","mea","j2r2","c22r2","spaop","spinc","splan","sprate"]
    
    columns=[]
    
    for objectnum in objects:
        for param in objparams:
            columns.append(param+"_"+objectnum)

    mmparamdf=pandas.DataFrame(index=[0],columns=columns)
    
    # values from ORBITS AND MASSES OF THE SATELLITES OF THE DWARF PLANET HAUMEA = 2003 EL61
    # Ragozzine & Brown (2009)
    
    names = ["Sun","Haumea","Namaka","Hi'iaka"]
    masses = [1.988e30, 4.006e21, 1.79e18, 1.79e19] # masses of the bodies in kg

    for objectnum in objects:
        mmparamdf["name_"+objectnum] = name
        mmparamdf["mass_"+objectnum] = masses[objectnum]
    
    # orbital parameters
    mmparamdf["sma0"] = 6.46e9  # km
    mmparamdf["ecc0"] = 0.194
    mmparamdf["aop0"] = np.pi/180.0*238.778 # radians
    mmparamdf["inc0"] = np.pi/180.0*28.214
    mmparamdf["lan0"] = np.pi/180.0*122.163
    mmparamdf["mea0"] = np.pi/180.0*217.774
    
    mmparamdf["sma2"] = 25657.0
    mmparamdf["ecc2"] = 0.249
    mmparamdf["aop2"] = np.pi/180.0*178.9
    mmparamdf["inc2"] = np.pi/180.0*13.013
    mmparamdf["lan2"] = np.pi/180.0*205.016
    mmparamdf["mea2"] = np.pi/180.0*178.5
    
    mmparamdf["sma3"] = 49880
    mmparamdf["ecc3"] = 0.0513
    mmparamdf["aop3"] = np.pi/180.0*154.1
    mmparamdf["inc3"] = np.pi/180.0*126.356
    mmparamdf["lan3"] = np.pi/180.0*206.766
    mmparamdf["mea3"] = np.pi/180.0*152.8
    
    mmparamdf["j2r2_1"] = 1.04e5 #km^2 
    mmparamdf["spaop_1"] = mmparamdf["aop3"] # Haumea's axis is aligned (roughly) with Hi'iaka's orbit
    mmparamdf["spinc_1"] = mmparamdf["inc3"]
    mmparamdf["splan_1"] = mmparamdf["lan3"]
    
    mmparamdf.fillna(0.0)
    
    return mmparamdf


#test_runprops()
print(test_mmparamdf())
