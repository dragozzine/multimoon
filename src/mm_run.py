# mm_run.py
# 
# The main script that runs MultiMoon
# Darin Ragozzine
# March 27, 2020
#
# Run with 
# python mm_run.py args
# where args are some combination of input and output file names TBD



# - Input: filename that is a JSON text file version of run_props dictionary
# - Output: directly, information about the success of the run
# - Output: indirectly: posterior, plots, etc.

#Works with run_props dictionary
# - Reads in run_props dictionary (from JSON)
# - Adds information autogenerated for this specific run
# - Checks other information

#Generates Initial Guess
#Make param_to_fit_scale: a parameters dataframe with the scale we'll use for converting to the fit values
#default: param_to_fit_scale is the first row of the init_guess dataframe
# - that way, every system and every parameter is pretty naturally scaled automatically!

#Checks if geocentric_object_position.csv already exists; if not, creates it; then reads it in and stores it #somewhere global

#Opens obsdata file and puts it into dataframe and stores it somewhere global

#Initializes Run
#Chooses prior and likelihood methods

#Runs emcee

#Takes posterior from emcee, thins, converts from fit_array to posterior dataframe in parameter format
#Saves posterior dataframe to designated output file

#Makes requested plots


"""Run MultiMoon

Inputs:
	filename that is a JSON text file version of run_props dictionary

Outputs:
	diagnostic information about the run

"""

import sys
import numpy as np
import pandas as pd
import emcee
import random
import h5py
#from tqdm import tqdm  # progress bar for emcee, but needs package
import mm_runprops
import mm_init_guess
import mm_likelihood
import mm_make_geo_pos
import mm_priors
import mm_relast
import mm_autorun
import mm_param
import os


# Read in the run props dictionary
# Adds information autogenerated for this specific run
# Checks other information
print('All import statements ran through.')

runprops = mm_runprops.runprops

nwalkers = runprops.get("nwalkers")

# Generate the intial guess for emcee
# starting guess is given by the user as specified in runprops
# and turned into the official initial guess

# LATER TODO: starting -> initial guess function is specificed by user

guesses = mm_init_guess.mm_init_guess(runprops)	# maybe more args
# ouptut from init_guess is a dataframe with all the desired parameters to be fit

# Getting relevant checking flags from runprops
dynamicstoincludeflags = runprops.get("dynamicstoincludeflags")
includesun = runprops.get("includesun")
paramnames = list(sum(list(guesses), ()))

# Check to make sure that numobjects equals length of dynamics flag
if len(dynamicstoincludeflags) != runprops.get("numobjects"):
    print("ERROR: Number of objects given in runprops.txt does not match the length of dynamicstoincludeflags")
    sys.exit()

# Now checking each object sequentially
for i in range(runprops.get("numobjects")):
    if i == 0:
        if dynamicstoincludeflags[i] == "0":
            if not (("mass_" + str(i+1) in paramnames)):
                print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                sys.exit()
            if (("sprate_" + str(i+1) in paramnames) or ("j2r2_" + str(i+1) in paramnames) or
		("spinc_" + str(i+1) in paramnames) or ("splan_" + str(i+1) in paramnames) or
		("c22r2_" + str(i+1) in paramnames) or ("spaop_" + str(i+1) in paramnames)):
                print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                sys.exit()
        elif dynamicstoincludeflags[i] == "1":
            if not (("mass_" + str(i+1) in paramnames) and ("j2r2_" + str(i+1) in paramnames) and
		    ("spinc_" + str(i+1) in paramnames) and ("splan_" + str(i+1) in paramnames)):
                print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                sys.exit()
            if (("sprate_" + str(i+1) in paramnames) or
		("c22r2_" + str(i+1) in paramnames) or ("spaop_" + str(i+1) in paramnames)):
                print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                sys.exit()
        elif dynamicstoincludeflags[i] == "2":
            if not (("mass_" + str(i+1) in paramnames) and ("j2r2_" + str(i+1) in paramnames) and
		    ("spinc_" + str(i+1) in paramnames) and ("splan_" + str(i+1) in paramnames) and
		    ("c22r2_" + str(i+1) in paramnames) and ("spaop_" + str(i+1) in paramnames) and
		    ("sprate_" + str(i+1) in paramnames)):
                print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                sys.exit()
        else:
            print("ERROR: dynamicstoincludeflags contains unallowed numbers. Allowed numbers are 0, 1, 2.")
            sys.exit()
    else:
        if dynamicstoincludeflags[i] == "0":
            if not (("mass_" + str(i+1) in paramnames) and ("sma_" + str(i+1) in paramnames) and
		    ("ecc_" + str(i+1) in paramnames) and ("inc_" + str(i+1) in paramnames) and
		    ("aop_" + str(i+1) in paramnames) and ("lan_" + str(i+1) in paramnames) and
		    ("mea_" + str(i+1) in paramnames)):
                print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                sys.exit()
            if (("sprate_" + str(i+1) in paramnames) or ("j2r2_" + str(i+1) in paramnames) or
		("spinc_" + str(i+1) in paramnames) or ("splan_" + str(i+1) in paramnames) or
		("c22r2_" + str(i+1) in paramnames) or ("spaop_" + str(i+1) in paramnames)):
                print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                sys.exit()
        elif dynamicstoincludeflags[i] == "1":
            if not (("mass_" + str(i+1) in paramnames) and ("sma_" + str(i+1) in paramnames) and
		    ("ecc_" + str(i+1) in paramnames) and ("inc_" + str(i+1) in paramnames) and
		    ("aop_" + str(i+1) in paramnames) and ("lan_" + str(i+1) in paramnames) and
		    ("mea_" + str(i+1) in paramnames) and ("j2r2_" + str(i+1) in paramnames) and
                    ("spinc_" + str(i+1) in paramnames) and ("splan_" + str(i+1) in paramnames)):
                print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                sys.exit()
            if (("sprate_" + str(i+1) in paramnames) or
		("c22r2_" + str(i+1) in paramnames) or ("spaop_" + str(i+1) in paramnames)):
                print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                sys.exit()
        elif dynamicstoincludeflags[i] == "2":
            if not (("mass_" + str(i+1) in paramnames) and ("sma_" + str(i+1) in paramnames) and
		    ("ecc_" + str(i+1) in paramnames) and ("inc_" + str(i+1) in paramnames) and
		    ("aop_" + str(i+1) in paramnames) and ("lan_" + str(i+1) in paramnames) and
		    ("mea_" + str(i+1) in paramnames) and ("j2r2_" + str(i+1) in paramnames) and
                    ("spinc_" + str(i+1) in paramnames) and ("splan_" + str(i+1) in paramnames) and
		    ("c22r2_" + str(i+1) in paramnames) and ("spaop_" + str(i+1) in paramnames) and
		    ("sprate_" + str(i+1) in paramnames)):
                print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                sys.exit()
        else:
            print("ERROR: dynamicstoincludeflags contains unallowed numbers. Allowed numbers are 0, 1, 2.")
            sys.exit()

# Now checking the includesun flag
if includesun:
    if not (("mass_0" in paramnames) and ("sma_0" in paramnames) and
	    ("ecc_0" in paramnames) and ("inc_0" in paramnames) and
	    ("aop_0" in paramnames) and ("lan_0" in paramnames) and
	    ("mea_0" in paramnames)):
        print("Error: includesun flag does not match inputs.")
        sys.exit()
if not includesun:
    if (("mass_0" in paramnames) or ("sma_0" in paramnames) or
	("ecc_0" in paramnames) or ("inc_0" in paramnames) or
	("aop_0" in paramnames) or ("lan_0" in paramnames) or
	("mea_0" in paramnames)):
        print("Error: includesun flag does not match inputs.")
        sys.exit()
    
#ndim is equal to the number of dimension, should this be equal to the number of columns of the init_guess array?
ndim = len(guesses.columns)

# Convert the guesses into fitting units and place in numpy array
p0,float_names,fixed_df,total_df_names,fit_scale = mm_param.from_param_df_to_fit_array(guesses,runprops)
#we still do not have a constraints or fit scale defined

# Check to see if geocentric_object_position.csv exists and if not creates it
objname = runprops.get('objectname')
if os.path.exists("../data/" + objname + "/geocentric_" + objname + "_position.csv"):
	print("Object geocentric position file geocentric_" + objname + "_position.csv will be used")
else:
	print("No object geocentric position file exists. Creating new file.")
	mm_make_geo_pos.mm_make_geo_pos(objname, start='2000-01-01', end='2040-01-01', step='10d')	# This is basically a function based on DS's makeHorFile
	print("geocentric_" + objname + "_position.csv has been created")

# Reads in th geocentric_object data file
geocentric_object_positions = pd.read_csv("../data/" + objname + "/geocentric_" + objname + "_position.csv")

# Now get observations data frame
# DS TODO: take observations data frame from runprops
obsdata = runprops.get('obsdata_file')

obsDF = 0
if os.path.exists(obsdata):
	print("Observational data file " + obsdata + " will be used")
	obsdf = pd.read_csv(obsdata)
else:
	print("ERROR: No observational data file exists. Aborting run.")
	sys.exit()

# Go through initial guesses and check that all walkers have finite posterior probability
reset = 0
maxreset = runprops.get("maxreset")
for i in range(nwalkers):  
	llhood = mm_likelihood.log_probability(p0[i,:], float_names,fixed_df,total_df_names, fit_scale, runprops, obsdf)
	while (llhood == -np.Inf):
		# Resetting walker to be average of two other walkers
		# BP TODO: Test this to make sure it works...
		p0[i,:] = (p0[random.randrange(nwalkers),:] + p0[random.randrange(nwalkers),:])/2
		llhood = mm_likelihood.log_probability(p0[i,:], float_names,fixed_df,total_df_names, fit_scale, runprops, obsdf)
		reset += 1
		if reset > maxreset:
			print("Maximum number of resets has been reached, aborting run.")
			sys.exit() 

# We now have an initial guess for each walker that is not really bad.
# Begin MCMC

# Now creating the sampler object
filename = "tutorial.h5"	# BP TODO: rename this to be something meaningful
# BP TODO: make an option in runprops to start from the end of another run and just append it
backend = emcee.backends.HDFBackend(filename)
backend.reset(nwalkers, ndim)

sampler = emcee.EnsembleSampler(nwalkers, ndim, 
	mm_likelihood.log_probability, backend=backend, 
	args = (runprops, fitarray_to_params_dict, obsdf))

#Starting the burnin
# BP TODO: autoburnin??
# So looking at how the emcee documentation does burn ins while saving the file, it seems like
# the best way to do a run is to just do a single long run until the ess > 100 and cut the 
# burn in off afterwards. This way you can save the whole chain and you can really analyze where to
# cut off the burn in.

nburnin = runprops.get("nburnin")
print("Starting the burn in")
state = sampler.run_mcmc(p0, nburnin, progress = True, store = False)
sampler.reset()

# Now do the full run with essgoal and initial n steps

nsteps = runprops.get("nsteps")
essgoal = runprops.get("essgoal")
maxiter = runprops.get("maxiter")

mm_autorun.mm_autorun(sampler, essgoal, state, initsteps, maxiter)

# Once it's completed, we need to save the chain
chain = sampler.get_chain(thin = runprops.get("nthinning"))
flatchain = sampler.get_chain(flat = True, thin = runprops.get("nthinning"))

# TODO: save chains to file

# make plots of MCMC results
mm_analysis.plots(sampler,guesses.columns)
mm_analysis.autocorrelation(sampler,guesses.columns)

# make other diagnostic plots
# TODO: orbit astrometry plots
# TODO: residual plots
