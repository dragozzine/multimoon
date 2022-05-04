#
#	mm_clustering.py
#	
#	Implements the clustering algorithm described in Hou 2010
#
#	Benjamin Proudfoot
#	08/04/20
#

# load in all required packages and functions
import numpy as np
import pandas as pd
import random
import commentjson as json
import emcee

def mm_clustering(sampler, state, float_names, fixed_df, total_df_names, fit_scale, runprops, obsdf,geo_obj_pos, best_llhoods, backend, pool, mm_likelihood, ndim, moveset, const = 50, lag = 10, max_prune_frac = 0.9):
	"""
	Determines if walkers in the ensemble are lost and removes them. Replacing them
	with random linear combinations of two random good wwalkers. 
	Inputs:
		sampler: an emcee sampler object on which to run MultiMoon

		state: ndarray, a place from which to start the emcee run

		float_names-fit_scale: internal multimoon data types used in mm_likelihood

		runprops: internal runprops dictionary

		obsdf, geo_obj_pos: internal observational data types

		best_llhood-pool: more internal inputs

		ndim: number of dimensions of parameter vector

		moveset: the moveset you wwant emcee to use

		const: constant used in Hou 2010 clustering algorithm. Determines how tolerant the clustering is.

		lag: the number of steps to average the likelihood over.

		max_prune_frac: maximum fraction of wwalkers to remove.

	Outputs:
		emcee sampler object, state object
	"""

	# Get inputs from runprops
	nwalkers = runprops.get("nwalkers")
	reburnin = runprops.get("clustering_burnin")
	if reburnin == 0:
		return sampler, state
	verbose = runprops.get("verbose")
	nthinning = runprops.get("nthinning")

	# Getting important values from the chain
	lastparams = sampler.get_chain()[-1,:,:]
	ngens = sampler.get_chain().shape[0]
	
	if ngens < lag:
		if verbose:
			print("Chain too short for clustering algorithm, clustering not performed")
		return sampler, state
	avllhood = np.mean(sampler.get_log_prob()[-lag:,:], axis = 0)

	if verbose:
		print(np.sort(avllhood))
		print(sampler.acceptance_fraction)

	# Sorting the walkers by likelihood values
	llhoodparam = pd.DataFrame(columns = ['llhood'] + float_names)
	for i in range(nwalkers):
		llhoodparam.loc[i] = np.concatenate([np.array([avllhood[i]]),lastparams[i,:]])
	llhoodparam = llhoodparam.sort_values(by=['llhood'], ascending = False)
	llhoodparam = llhoodparam.values

	# Performing rejection tests
	reject = np.zeros(nwalkers)
	for i in range(1,nwalkers-1):
		term1 = -llhoodparam[i+1,0] + llhoodparam[i,0]
		term2 = const*(-llhoodparam[i,0] + llhoodparam[0,0])/(i)
		print(term1, term2)
		if term1 > term2:
			reject[(i+1):] = 1
			break
	freject = reject.sum()/nwalkers
	if verbose:
		print(freject)

	# Pruning walkers based on the clusters found, 
	# replacing them with random linear combinations of walkers within the cluster
	# Skipping if cluster is not big enough
	if freject < max_prune_frac:
		params = llhoodparam[:,1:]
		for i in range(len(reject)):
			if reject[i] == 1:
				# Picking random inputs
				p = random.random()
				c1 = random.randrange(i)
				c2 = random.randrange(i)
				while c1 == c2:
					c2 = random.randrange(i)
				# Calculating random linear comb of two random walkers
				params[i,:] = (p*params[c1,:] + (1-p)*params[c2,:])
		# Create a new sampler
		sampler = emcee.EnsembleSampler(nwalkers, ndim, mm_likelihood.log_probability, 
						backend=backend, pool = pool,
						args = (float_names, fixed_df, total_df_names, fit_scale, runprops, obsdf,geo_obj_pos, best_llhoods),
						moves = moveset)
		# Perform another burn in
		if runprops.get('thin_run'):
			state = sampler.run_mcmc(params, reburnin, progress = True, store = True, thin=nthinning)
		else:
			state = sampler.run_mcmc(params, reburnin, progress = True, store = True)
		return sampler, state
	else:
		print("Cluster not big enough, clustering not performed")
		return sampler, state
