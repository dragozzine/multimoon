#
#	scrape.py
#
#	Scrapes all of the completed Keplerian runs for the 2021 227 class project
#
#	Benjamin Proudfoot
#	10/13/21
#

import numpy as np
import pandas as pd
import os
from tqdm import tqdm
from scipy import interpolate
import sys
import commentjson as json

class ReadJson(object):
    def __init__(self, filename):
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

# What do we want in the output summary?
# -best log likelihood
# -best reduced chi squared
# -best p-value
# -Best fits
# -Parameter distributions
#	-Median
#	-Sigmas
#	-How do we show a detection of J2 in a table....
#		-Maybe we don't have to... students can identify a J2 detection by eye (e.g. Borasisi)
#
# -have this code also make ensemble figures???
#	-Overplot all J2 posteriors. How to deal with mirror orbits??
#	-Begin thinking about spin-orbit angle plots... 
#		-How do you display non-detections (uniform distribution) with detections.
#
#
#
# Load in the names of the results directories
results = np.loadtxt("resultfiles_locked_2.txt", delimiter = "\n", dtype = "str")
rows = results.copy()
for i,file in enumerate(results):
	results[i] = file[6:]
	rows[i] = file[17:]

# Loop over results and find zscores
table = pd.DataFrame(index = rows, columns = ["lhood","reduced_chi_sq","mass_tot","sma","ecc","aop","inc","lan","mea","j2r2","mass_1_best","mass_2_best","sma_best","ecc_best","aop_best","inc_best","lan_best","mea_best","j2r2_best"])
table = table.transpose().copy()

for i,file in enumerate(results):
	if not os.path.isdir(file):
		print(file + " doesn't exist yet :(")
		continue

	print(file)

	runprops = ReadJson(file + "runprops.txt").outProps()
	objname = runprops.get("objectname")

	# Getting values from the best_likelihoods file
	bestlfile = pd.read_csv(file + "best_likelihoods.csv", index_col = False)
	#print(bestlfile.iloc[-1,:])
	bestl = bestlfile.iloc[-1,0]
	bestrcs = bestlfile.iloc[-1,4]
	dof = bestlfile.iloc[-1,1]
	pval = bestlfile.iloc[-1,2]
	bestparams = bestlfile.iloc[-1,6:-5]
	#print(bestparams)

	# Getting Grundy best fits
	truevals = pd.read_csv(file + objname + "_init_guess.csv", sep = ",")
	truths = []
	truths.append(truevals.values[0,1] + truevals.values[1,1])
	truths.append(truevals.values[2,1])
	truths.append(truevals.values[3,1])
	truths.append(truevals.values[4,1])
	truths.append(truevals.values[5,1])
	truths.append(truevals.values[6,1])
	truths.append(truevals.values[7,1])
	truths.append(truevals.values[8,1])
	truths = np.array(truths).astype(np.float)

	# Getting sigsdf
	try:
		sigsdf = pd.read_csv(file + "sigsdf.csv", sep = ",").values
	except:
		print(file + "sigsdf.csv doesn't exist yet :(")
		continue

	sigsdf = sigsdf[:,1:-2]
	for k in range(7):
		sigsdf[k,:] = sigsdf[k,:] + sigsdf[k,3]
	sigsdf[:,3] = sigsdf[:,3]/2
	sigsdf[1,:] = sigsdf[0,:] + sigsdf[1,:]
        sigsdf = copy(sigsdf[1:,:])
	print(sigsdf)
	print()

	#Calculating zscores
	zscores = []
	sigma = np.array([-3,-2,-1,0,1,2,3]).astype(np.float)
	for k in range(7):
		vals = sigsdf[k,:].flatten().astype(np.float)
		#print(vals)
		f = interpolate.interp1d(vals,sigma, kind = "cubic", bounds_error = False,
					 fill_value = "extrapolate")
		z = f(truths[k])
		#print(vals, truths[k])
		if z > 3.0:
			z = ">3"
		elif z < -3.0:
			z = "<-3"
		zscores.append(z)
	row = np.array([bestl])
	row = np.append(row, [bestrcs])
	row = np.append(row, [pval])
	row = np.append(row, zscores)
	row = np.append(row, bestparams)
	table.iloc[:,i] = row

table = table.transpose().copy()
table.to_csv("../results/scraped_results_locked.csv")
