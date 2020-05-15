#
#	mm_autorun.py
#
#	The purpose of this program is to autorun multimoon until an ESS (effective sample
#	size) is reached.
#
#	Benjamin Proudfoot
#	05/15/20
#

import sys
import numpy as np
import emcee
from tqdm import tqdm

def mm_autorun(essgoal, sampler, state, initsteps, maxiter):
	
	# First run the sampler for the specified number of initial steps
	sampler.run_mcmc(state, initsteps)
