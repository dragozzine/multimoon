#
#	scrape_posteriors.py
#
#	Scrapes all chain.h5 files to create plots
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
import emcee
from unpack import *
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
results = np.loadtxt("resultfiles_fulllong.txt", delimiter = "\n", dtype = "str")
rows = results.copy()
for i,file in enumerate(results):
	results[i] = file[6:]
	rows[i] = file[17:]

# Loop over results and create aggregated plots

# Create plots be editing
j2fig = plt.figure()
angfig = plt.figure()

for i,file in enumerate(results):
	if not os.path.isdir(file + "sigsdf.csv"):
		print(file + " hasn't finished yet :(")
		continue

	print(file)
	runprops = ReadJson(file + "runprops.txt").outProps()

	# Unpack chain file
	flatchain, names = unpack(file, thin_plots = 100)

	# Isolate relevant indices in the flatchain
	j2r2_index = [n for n, l in enumerate(names) if l.startswith('j2r2_1')][0]
	spinc1_index = [n for n, l in enumerate(names) if l.startswith('spinc_1')][0]
	splan1_index = [n for n, l in enumerate(names) if l.startswith('splan_1')][0]
	inc_index = [n for n, l in enumerate(names) if l.startswith('inc_2')][0]
	lan_index = [n for n, l in enumerate(names) if l.startswith('lan_2')][0]
	
	j2_arr = np.log10(np.exp(flatchain[:,j2r2_index])/(runprops.get("c_axis")**2))
	spinc1_arr = np.deg2rad(flatchain[:,spinc1_index])
	splan1_arr = np.deg2rad(flatchain[:,splan1_index])
	inc_arr = np.deg2rad(flatchain[:,inc_index])
	lan_arr = np.deg2rad(flatchain[:,lan_index])

	# Adding J2 histogram
	j2fig.hist(j2_arr, density = True, histtype = "step")

	# Calculating and adding spin orbit angle
	mutualinc = np.rad2deg(np.arccos( np.cos(spinc1_arr)*np.cos(inc_arr) + np.sin(spinc1_arr)*np.sin(inc_arr)*np.cos(splan1_arr - lan_arr) ))
	angfig.hist(mutualinc, density = True, histtype = "step")


j2fig.xlabel("log10(J2)")
j2fig.savefig("j2post.pdf")
angfig.xlabel("Spin-Orbit Mutual Inclination (deg)")
angfig.savefig("angpost.pdf")
