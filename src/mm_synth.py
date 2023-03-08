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

# load in all required packages and functions
import numpy as np
import pandas as pd
import sys
import h5py
import random
import commentjson as json
import os
import shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


import mm_priors as prior
import mm_likelihood
import mm_param
import mm_runprops
import mm_make_geo_pos

sys.path.insert(1, 'mm_SPINNY/')
from mm_SPINNY.spinny_vector import generate_vector
import mm_relast

# Start getting things needed
runprops = mm_runprops.runprops

verbose = runprops.get("verbose")
nobjects = runprops.get("numobjects")
objname = runprops.get('objectname')

# Setting the observations data file and geo position data file
obsdata = runprops.get('obsdata_file')

obsdf = 0
if os.path.exists(obsdata):
	if verbose:
		print("Observational data file " + obsdata + " will be used")
	obsdf = pd.read_csv(obsdata)
else:
	print("ERROR: No observational data file exists. Aborting run.")
	sys.exit()

times = obsdf['time'].tolist()
geo_obj_pos = mm_make_geo_pos.mm_make_geo_pos(objname, times, runprops, True)

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

#getData = ReadJson("../runs/runprops.txt")
#runprops = getData.outProps()

# Outputting model astrometry based on the params df
Model_DeltaLong, Model_DeltaLat, obsdf = mm_likelihood.mm_chisquare(paramdf, obsdf, runprops, geo_obj_pos, gensynth = True)

# Changing the astrometry based on the effects of a hidden moon.
if runprops.get('photo_offset'):
        
    mass_ratio = paramdf['mass_2'][0]/paramdf['mass_1'][0]
     
    f_val = paramdf['f_val_1'][0]
    
    bright_ratio = f_val*mass_ratio**(2/3)
        
    rel_pos_lat = Model_DeltaLat[0,:]
    rel_pos_long = Model_DeltaLong[0,:]
    print("Models: ",Model_DeltaLat)
    print(Model_DeltaLong)

    delta_offset_lat = bright_ratio*rel_pos_lat
    delta_offset_long = bright_ratio*rel_pos_long
        
    Model_DeltaLat = Model_DeltaLat + delta_offset_lat
    Model_DeltaLong = Model_DeltaLong + delta_offset_long

# Putting the synthetic data into the right formats
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

# Adding random errors to obsdf
q = 1.0
for i in range(1,nobjects):
	obsdf["DeltaLat_" + objectnames[i]] = obsdf["DeltaLat_" + objectnames[i]].values + q*np.random.normal(size = obsdf.shape[0])*obsdf["DeltaLat_" + objectnames[i] + "_err"].values
	obsdf["DeltaLong_" + objectnames[i]] = obsdf["DeltaLong_" + objectnames[i]].values + q*np.random.normal(size = obsdf.shape[0])*obsdf["DeltaLong_" + objectnames[i] + "_err"].values

# Adding additional random errors to the data to test robust statistics
true_frac = 1.0
bkg_std = 0.01       # in arcseconds
bkg_lon = np.random.rand(obsdf.shape[0]) > true_frac
bkg_lat = np.random.rand(obsdf.shape[0]) > true_frac

for i in range(1,nobjects):
	obsdf["DeltaLat_" + objectnames[i]][bkg_lat] += bkg_std*np.random.randn(sum(bkg_lat))
	obsdf["DeltaLong_" + objectnames[i]][bkg_lon] += bkg_std*np.random.randn(sum(bkg_lon))

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
savedir = "../runs/" + runprops.get("run_dir") + "/" + runprops.get("runs_file") + "/"

try:
	os.mkdir(savedir)
	print("Made directory to save testcase (" + savedir + ").")
except FileExistsError:
	print("directory already exists. Removing obsdf and .png files")
plt.savefig(savedir + "astrometry.png")
obsdf.to_csv(savedir + runprops.get("run_dir")+"_obs_df.csv")
shutil.copy("runprops_gensynth.txt", savedir)
geo_obj_pos.to_csv(savedir + runprops.get("run_dir")+"geocentric_pos.csv")
