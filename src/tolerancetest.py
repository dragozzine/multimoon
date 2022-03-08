#
#	tolerancetest.py
#
#	Tests to see what SPINNY tolerance should be used for your models
#	
#	Benjamin Proudfoot
#	03/02/2021
#
#	To use, replace the init files below with your desired run files. In addition, input the
#	intervals you want tolerance tested at in 'tolvals'. The script will take these and 
#	create a single "correct" model with tolerance=1e-15, and then create subsequent models
#	at the tolvals. These models will be compared to the reference "correct" model and
#	a delta chi-squared will be calculated. The comparisons, along with times to run each model,
#	will be output into two files. We recommend using a tolerance which keeps the delta 
#	chi-squared below 1, while maintaining the shortest possible run time. 
#

objectname = "2006 BR284"
runtype = "20"

# Code starts
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
#import mm_runprops
import mm_make_geo_pos
sys.path.insert(1, 'mm_SPINNY/')
from mm_SPINNY.spinny_vector import generate_vector
import mm_relast
import time

# Input files
print("../runs/" + objectname + "/" + runtype + "/runprops.txt")
runprops = ReadJson("../runs/" + objectname + "/" + runtype + "/runprops.txt").outProps()
initparams = pd.read_csv("../runs/" + objectname + "/" + runtype + "/" + objectname + "_init_guess.csv", index_col = 0)
#obsdata = runprops.get('
obsdf = runprops.get('obs_df')
obsdata = "../runs/" + objectname + "/observations/" + obsdf
#runprops = ReadJson("../data/227_2021/runs/2006 BR284/10_locked/runprops.txt").outProps()
#initparams = pd.read_csv("../data/227_2021/runs/2006 BR284/10_locked/2006 BR284_init_guess.csv", index_col = 0)
#obsdata = "../data/227_2021/runs/2006 BR284/observations/2006 BR284_obs_df.csv"

tolvals = np.logspace(-14,-8,20)

verbose = runprops.get("verbose")
nobjects = runprops.get("numobjects")
objname = runprops.get('objectname')

params = []
paramnames = []
params2 = []
paramnames2 = []
names = []

name_dict = runprops.get("names_dict")

for k in initparams.iloc[:,0].values:
	params.append(k)
	params2.append(k)
for k in list(initparams.index):
	paramnames.append(k)
	paramnames2.append(k)
if runprops.get('includesun') == 1:
	paramnames.append('name_0')
	paramnames2.append('name_0') 
#params.append(500.0)
#params2.append(500.0)
#params.append(1.0)
#params2.append(1.0)
#paramnames.append("ax_1")
#paramnames2.append("ax_1")
if runprops.get('includesun') == 1:
	params.append('Sun')
	params2.append('Sun')    
for k in name_dict.values():
	params.append(k)
	params2.append(k)
	names.append(k)
for k in name_dict.keys():
	paramnames.append(k)
	paramnames2.append(k)

paramdf = pd.DataFrame(params).transpose()
paramdf.columns = paramnames
paramdf2 = pd.DataFrame(params2).transpose()
paramdf2.columns = paramnames2

print(paramdf.to_string())
print(paramdf2.to_string())

obsdf = 0
if os.path.exists(obsdata):
	if verbose:
		print("Observational data file " + obsdata + " will be used")
	obsdf = pd.read_csv(obsdata)
else:
	print(obsdata)
	print("ERROR: No observational data file exists. Aborting run.")
	sys.exit()

times = obsdf['time'].tolist()
geo_obj_pos = mm_make_geo_pos.mm_make_geo_pos(objname, times, runprops, True)

#print(times)
start_time = time.time()
runprops["spinny_tolerance"] = 1e-15
#print('paramdf',paramdf)
Model_DeltaLong, Model_DeltaLat, obsdf = mm_likelihood.mm_chisquare(paramdf, obsdf, runprops, geo_obj_pos, gensynth = True)
print("--- %s seconds ---" % (time.time() - start_time))

# Loop over tolerance values
dchisq = np.zeros(tolvals.size)
runtimes = np.zeros(tolvals.size)
for k, tol in enumerate(tolvals):
	start_time = time.time()
	runprops["spinny_tolerance"] = tol
	Model_DeltaLong2, Model_DeltaLat2, obsdf = mm_likelihood.mm_chisquare(paramdf2, obsdf, runprops, geo_obj_pos, gensynth = True)
	t = time.time() - start_time
	print("Tolerance = ", tol)
	print("--- %s seconds ---" % (t))
	runtimes[k] = t

	rows = obsdf.shape[0]
	residuals = np.zeros(((nobjects-1)*2, rows))

	for i in range(rows):
		for j in range(1,nobjects):
			# Add noise
			#Model_DeltaLong[j-1][i] = Model_DeltaLong[j-1][i] + np.random.normal()*obsdf["DeltaLong_" + names[j] + "_err"][i]
			#Model_DeltaLat[j-1][i] = Model_DeltaLat[j-1][i] + np.random.normal()*obsdf["DeltaLat_" + names[j] + "_err"][i]
	
			residuals[2*(j-1)][i] = ((Model_DeltaLong[j-1][i]-Model_DeltaLong2[j-1][i])/obsdf["DeltaLong_"+names[j]+"_err"][i])
			residuals[2*(j-1)+1][i] = ((Model_DeltaLat[j-1][i]-Model_DeltaLat2[j-1][i])/obsdf["DeltaLat_"+names[j]+"_err"][i])
			#print(names[j])
			#print(residuals[2*(j-1)][i],residuals[2*(j-1)+1][i])

	chisquares = residuals**2
	chisq_tot = np.zeros(2*nobjects)
	for i in range(0,2*nobjects-2):
		chisq_tot[i]=np.nansum(chisquares[i])
	chisquare_total = np.nansum(chisq_tot)
	print(residuals)
	print(chisquare_total)	
	dchisq[k] = chisquare_total

# Plotting "error values"
plt.figure()
plt.loglog(tolvals, dchisq, marker = ".")
plt.xlabel("SPINNY Tolerance Value")
plt.ylabel(r"$\Delta \chi^2$ (compared to tol=1e-15")
plt.savefig("../runs/" + objectname + "/" + runtype + "/toltest.png")
plt.close()

plt.figure()
plt.loglog(tolvals, runtimes, marker = ".")
plt.xlabel("SPINNY Tolerance Value")
plt.ylabel("Runtime (s)")
plt.savefig("../runs/" + objectname + "/" + runtype + "/toltest_time.png")
plt.close()



