#
#	mm_autorun.py
#
#	The purpose of this program is to autorun multimoon until an ESS (effective sample
#	size) is reached.
#
#	Benjamin Proudfoot
#	05/15/20
#
#	TODO: Add some methods which can estimate the time?? or maybe just use tqdm
#
#
#
#
#

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import numpy as np
import emcee
from tqdm import tqdm
from mm_analysis import *

def mm_autorun(sampler, essgoal, state, initsteps, maxiter):
	"""
	Runs MultiMoon until an effective sample size goal is achieved
	or the maximum number of iterations is reached. 

	Inputs:
		sampler: an emcee sampler object on which to run MultiMoon

		essgoal: int, an effective sample size that the user wishes to reach

		state: ndarray, a place from which to start the emcee run

		initsteps: int, an initial guess as to how many steps are needed

		maxiter: int, a maximum number of steps to take

	Outputs:
		emcee sampler object

	"""
	# First run the sampler for the specified number of initial steps
	print("Running MultiMoon for ", initsteps, " steps")
	state = sampler.run_mcmc(state, initsteps, progress = True)

	# Find the integrated autocorrelation time (iat) of the run
	print(initsteps, " steps have been completed.")
	print("Calculating the effective sample size.")
	iat = autocorrelation(sampler, filename = "_0")

	# Assess the accuracy of the iat estimate
	ngens = sampler.get_chain().shape[0]	
	counter = 1
	while (50*iat > ngens):
		print("Runnning the sampler for another ", initsteps, " steps")
		state = sampler.run_mcmc(state, initsteps, progress = True)
		print(initsteps, " steps have been completed.")
		print("Calculating the effective sample size.")
		iat = autocorrelation(sampler, filename = "_" + str(counter))
		ngens = sampler.get_chain().shape[0]	
		counter += 1
		if ngens > maxiter:
			print("Maximum iterations has been reached, ending automated runs.")
			print("WARNING: The estimate for ESS cannot be trusted.")
			return sampler, ngens/iat

	# Now we have a reliable estimate for the iat, moving on to find ess		
	ess = ngens/iat
	
	# Ending if we have achieved goal
	if ess > essgoal:
		return sampler, ess

	# Calculating the number of steps to reach the ess goal
	nadditional = int(iat*essgoal - ngens)

	# Adding in a 10% buffer on nadditional to be sure goal is reached
	nadditional = int(nadditional*1.1)

	# Making sure maxiter is not reached
	if int(nadditional + ngens) > maxiter:
		nadditional = int(maxiter - ngens)
		state = sampler.run_mcmc(state, nadditional, progress = True)
		iat = autocorrelation(sampler, filename = "_" + str(counter))
		counter += 1
		ngens = sampler.get_chain().shape[0]
		print("Maximum iterations has been reached, ending automated runs.")
		return sampler, ngens/iat
	else:
		state = sampler.run_mcmc(state, nadditional, progress = True)
		iat = autocorrelation(sampler, filename = "_" + str(counter))
		counter += 1
		ngens = sampler.get_chain().shape[0]
		print("ESS goal has been reached, ending automated runs.")
		return sampler, ngens/iat




