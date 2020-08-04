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
import os
import shutil
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt


import mm_priors as prior
import mm_likelihood
import mm_param
import mm_runprops

sys.path.insert(1, 'mm_SPINNY/')
from mm_SPINNY.spinny_vector import generate_vector
import mm_relast

plt.switch_backend("MacOSX")

# maybe have a runprops type thing to hold the inputs?
# it then saves a copy of this to the directory where the test case is stored

# Getting the inputs from runprops_gensynth.txt
# This file needs to have everything that mm_chisquare needs
#	numobjects, verbose, float_dict, 
#filename = "runprops_gensynth.txt"

#getData= ReadJson(filename)
#runprops = getData.outProps()

runprops = mm_runprops.runprops


verbose = runprops.get("verbose")
nobjects = runprops.get("numobjects")

# Setting the observations data file and geo position data file
runprops["obsdata_file"] = "../data/" + runprops.get("objectname") + "/" + runprops.get("objectname") + "_obs_df.csv"
obsdata = runprops.get('obsdata_file')

obsdf = 0
if os.path.exists(obsdata):
	if verbose:
		print("Observational data file " + obsdata + " will be used")
	obsdf = pd.read_csv(obsdata)
else:
	print("ERROR: No observational data file exists. Aborting run.")
	sys.exit()


objname = runprops.get('objectname')
if os.path.exists("../data/" + objname + "/geocentric_" + objname + "_position.csv"):
	if verbose:
		print("Object geocentric position file geocentric_" + objname + "_position.csv will be used")
else:
	if verbose:
		print("No object geocentric position file exists. Aborting Run.")
	sys.exit()
geo_obj_pos = pd.read_csv("../data/" + objname + "/geocentric_" + objname + "_position.csv")

# Package the parameters wanted into a guesses-like df
params = []
paramnames = []
objectnames = []

params_dict = runprops.get("params_dict")
name_dict = runprops.get("names_dict")

for i in params_dict.values():
	params.append(i)
for i in params_dict.keys():
	paramnames.append(i)
for i in name_dict.values():
	params.append(i)
	objectnames.append(i)
for i in name_dict.keys():
	paramnames.append(i)
paramdf = pd.DataFrame(params).transpose()
paramdf.columns = paramnames

# Outputting model astrometry based on the params df
Model_DeltaLong, Model_DeltaLat, obsdf = mm_likelihood.mm_chisquare(paramdf, obsdf, runprops, geo_obj_pos, gensynth = True)
#positionData = mm_likelihood.mm_chisquare(paramdf, obsdf, runprops, geo_obj_pos, gensynth = True)

print(obsdf["time"])

cols = ["time","Lat_Prim","Long_Prim"]
for i in range(1,nobjects):
	obsdf["DeltaLat_" + objectnames[i]] = Model_DeltaLat[i-1]
	obsdf["DeltaLong_" + objectnames[i]] = Model_DeltaLong[i-1]
	cols.append("DeltaLat_" + objectnames[i])
	cols.append("DeltaLong_" + objectnames[i])
	cols.append("DeltaLat_" + objectnames[i] + "_err")
	cols.append("DeltaLong_" + objectnames[i] + "_err")

for i in range(obsdf.shape[0]):
	row = obsdf.iloc[i,:]
	for j in range(1,nobjects):
		if np.isnan(row["DeltaLat_" + objectnames[j] + "_err"]):
			obsdf.iloc[i,:]["DeltaLat_" + objectnames[j]] = np.nan
		if np.isnan(row["DeltaLong_" + objectnames[j] + "_err"]):
			obsdf.iloc[i,:]["DeltaLong_" + objectnames[j]] = np.nan

obsdf = obsdf.drop(labels = [col for col in obsdf if col not in cols], axis = 1)

# Now plot it to check to see if it look okay
x = np.empty((nobjects-1, obsdf.shape[0]))
xe = np.empty((nobjects-1, obsdf.shape[0]))
y = np.empty((nobjects-1, obsdf.shape[0]))
ye = np.empty((nobjects-1, obsdf.shape[0]))

fmts = ["bo","ro"]

fig = plt.figure()
for i in range(1,nobjects):
	x[i-1,:] = obsdf["DeltaLat_" + objectnames[i]].values
	xe[i-1,:] = obsdf["DeltaLat_" + objectnames[i] + "_err"].values
	y[i-1,:] = obsdf["DeltaLong_" + objectnames[i]].values
	ye[i-1,:] = obsdf["DeltaLong_" + objectnames[i] + "_err"].values
	plt.errorbar(x[i-1,:], y[i-1,:], xerr = xe[i-1,:], yerr = ye[i-1,:], fmt = fmts[i-1])

plt.axis('equal')

# Saving things now
savename = runprops.get("testcase_name")
savedir = "../testcases/" + savename + "/"

try:
	os.mkdir(savedir)
	if verbose:
		print("Made directory to save testcase (" + savedir + ").")
except FileExistsError:
	if verbose:
		print("directory already exists. Removing obsdf and .png files")
	os.remove(savedir + "astrometry.png")
	os.remove(savedir + "obsdf.csv")
plt.savefig(savedir + "astrometry.png")
obsdf.to_csv(savedir + "obsdf.csv")
shutil.copy("runprops_gensynth.txt", savedir)


# Need to think about excluding data points where the secondary/tertiary/etc is aligned with the 
# primary. In theory the fitter will have no problem fitting to this data, but it really isn't 
# realistic...

# Once the model dfs are output, then all that needs to be done is to output it as a data file.

# I think the best way to do all of this is to use a REAL observations file to get the date/time
# of observations and the error bars, and then replace the data points with synthetic ones.
