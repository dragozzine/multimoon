#
#	mm_synth.py
#
#	Generates synthetic astrometry to test multimoon
#
#	Benjamin Proudfoot
#	06/19/20
#

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

# I think the best method would be to have a tests directory where we store the synthetic astrometry,
# an output file associated with the creation of that synthetic astrometry (basically the parameters),
# and everything else needed to test the system (geocentric_obj_position.csv, etc).

# Starting psuedocode

# load in all required packages and functions
import numpy as np
import pandas as pd
import sys
import h5py
import random
import json

import mm_priors as prior
import mm_likelihood
import mm_param

sys.path.insert(1, 'mm_SPINNY/')
from mm_SPINNY.spinny_vector import generate_vector
import mm_relast

# maybe have a runprops type thing to hold the inputs?
# it then saves a copy of this to the directory where the test case is stored

# Getting the inputs from runprops_gensynth.txt
# This file needs to have everything that mm_chisquare needs
#	numobjects, verbose, float_dict, 
filename = "runprops_gensynth.txt"

getData= ReadJson(filename)
runprops = getData.outProps()

# Setting the observations data file
runprops["obsdata_file"] = # TODO: insert observations file here

# Package the parameters wanted into a guesses-like df

# Creating a paramdf
float_params,float_names,fixed_df,total_df_names,fit_scale = mm_param.from_param_df_to_fit_array(guesses,runprops)
name_dict = runprops.get("names_dict")
paramdf = mm_param.from_fit_array_to_param_df(float_params, float_names, fixed_df, total_df_names, fit_scale, name_dict)

# Outputting model astrometry based on the params df
Model_DeltaLong, Model_DeltaLat = mm_chisquare(paramdf, obdsf, runprops, gensynth = True)

# Need to think about excluding data points where the secondary/tertiary/etc is aligned with the 
# primary. In theory the fitter will have no problem fitting to this data, but it really isn't 
# realistic...

# Once the model dfs are output, then all that needs to be done is to output it as a data file.

# I think the best way to do all of this is to use a REAL observations file to get the date/time
# of observations and the error bars, and then replace the data points with synthetic ones.
