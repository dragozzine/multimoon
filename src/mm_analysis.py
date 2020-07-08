#	mm_analysis.py
#
#	Runs all of the analysis on the chains from mm_run
#
#	Benjamin Proudfoot
#	05/05/20
#

# To be honest I don't think that mm_run should make the chain a df. We would convert it to one
# and then immediately convert it back for the corner plots and walker plots. I don't think 
# the step is necessary esp with the huge size of the chain.

# WARNING: The indices for all of the work with the chains might be reversed to what I have here
# I used the same as those in all of my julia codes, and I don't know if the indices are reversed
# with emcee compared to AffineInvariantMCMC in julia. SO beware. 



import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import corner
import numpy as np
import emcee
import sys
from mm_runprops import runprops


#chain = (nwalkers, nlink, ndim)
def plots(sampler, parameters, objname, fit_scale, float_names):
			# Here parameters is whatever file/object will have the run params
	flatchain = sampler.get_chain(flat = True)
	fit = []

	for i in fit_scale.columns:
		if i[0] in float_names:
			val = fit_scale[i][0]
			fit.append(val)

#First fit the flatchain with the fit parameters    
	fchain = [[0] * len(flatchain[0])] * len(flatchain)    
	for i in range(len(flatchain)):
		row = []
		for j in range(len(flatchain[0])):
			val = flatchain[i][j]*fit[j]
			row.append(val)
		fchain[i] = row

	flatchain = np.array(fchain)

	chain = sampler.get_chain(flat = False)

#Now fit the chain 
	cchain = [[[0] * len(chain[0][0])] * len(chain[0])] * len(chain)    
	for i in range(len(cchain)):
		for j in range(len(cchain[0])):
			row = []
			for k in range(len(cchain[0][0])):
				val = chain[i][j][k]*fit[k]
				cchain[i][j][k] = val
#				print(val)
			print('\ncchain ',cchain[i][j],'\nchain ', chain[i][j])
#			cchain[i][j] = row


	cchain = np.array(cchain)
	print('chain: ', len(chain), len(chain[0]), len(chain[0][0]), len(cchain), len(cchain[0]), len(cchain[0][0]))
#	print(chain[:,0,0], '\n',cchain[:,0,0],'\n',chain[:,1,0], '\n',cchain[:,1,0])
#	print(chain)
	names = []
	for i in float_names:
		names.append(i)

	fig = corner.corner(flatchain, bins = 40, labels = names, show_titles = True, 
			    plot_datapoints = False, color = "blue", fill_contours = True,
			    title_fmt = ".4e")
	fig.tight_layout(pad = 1.08, h_pad = 0, w_pad = 0)
	for ax in fig.get_axes():
		ax.tick_params(axis = "both", labelsize = 20, pad = 0.5)
	fname = "../runs/"+objname+"_"+runprops.get("date")+"/corner.pdf"       
	fig.savefig(fname, format = 'pdf')
	plt.close("all")
	
	# Now make the walker plots
	numparams = chain.shape[2]
	numwalkers = chain.shape[1]
	numgens = chain.shape[0]
	for i in range(numparams):
		plt.figure()
		for j in range(numwalkers):
			plt.plot(np.reshape(chain[0:numgens,j,i], numgens))
		plt.ylabel(names[i])
		plt.xlabel("Generation")
		plt.savefig("../runs/"+objname+"_"+runprops.get("date")+"/walker_"+names[i]+".png")
		plt.close()

	# Likelihood plots
	llhoods = sampler.get_log_prob(flat = True)
	for i in range(numparams):
		plt.figure(figsize = (9,9))
		plt.subplot(221)
		for j in range(numwalkers):
			print('indices: ',i,j)
			print('chain: ',chain[:,j,i])
			print('cchain: ',cchain[:,j,i])
			plt.hist(chain[:,j,i], bins = 40, histtype = "step",
				color = "black",
				alpha = 0.4, density = True)
		plt.hist(chain[:,:,i].flatten(), bins = 40, histtype = "step", color = "black", density = True)
		plt.subplot(223)
		plt.scatter(flatchain[:,i].flatten(), llhoods.flatten(),
			    c = np.mod(np.linspace(0,llhoods.size - 1, llhoods.size), numwalkers),
			    cmap = "nipy_spectral", edgecolors = "none")
		plt.xlabel(names[i])
		plt.ylabel("Log(L)")
		plt.subplot(224)
		plt.hist(llhoods.flatten(), bins = 40, orientation = "horizontal", 
			 histtype = "step", color = "black")
		plt.savefig("../runs/"+objname+"_"+runprops.get("date")+"/likelihood.pdf", format = 'pdf')
		plt.close("all")

def auto_window(taus, c):
	m = np.arange(len(taus)) < c * taus
	if np.any(m):
		return np.argmin(m)
	return len(taus) - 1

def autocorr_new(y, c = 5.0):
	f = np.zeros(y.shape[1])
	for yy in y:
		f += emcee.autocorr.function_1d(yy)
	f /= len(y)
	taus = 2.0 * np.cumsum(f) - 1.0
	window = auto_window(taus, c)
	return taus[window]

def autocorrelation(sampler, objname, filename = "", thin = 1):
	# Getting chain for the first parameter to calculate different values
	chain = sampler.get_chain(thin = thin)
	
	nwalkers = sampler.nwalkers
	ndims = sampler.ndim
	nsteps = chain.shape[0]
	
	# Calculating values to calculate tau for
	# This chould be changed eventually
	N = np.exp(np.linspace(np.log(100), np.log(nsteps), 10)).astype(int)

	# Setting up array for tau estimates
	tau = np.empty( (len(N), ndims) )

	# Calculating tau for each value in N for each parameter
	for i in range(ndims):
		thischain = chain[:,:,i].T
		for j, n in enumerate(N):
			tau[j,i] = autocorr_new(thischain[:, :n])

	# Setting up to plot curves for tau in a grid
	x = 3
	y = ndims
	nrows = 1
	ncols = 3
	while x < y:
		y = y - x
		nrows += 1

	if ncols > ndims:
		ncols = ndims
	return np.mean(sampler.get_autocorr_time(quiet = True))
'''   
	fig, ax = plt.subplots(nrows = nrows, ncols = ncols, sharey=True, 
			       gridspec_kw={'wspace': 0},
			       figsize = (6.4*(ncols),4.8*(nrows)), 
			       squeeze = False)
	fig.suptitle("Autocorrelation estimates")
	fig.add_subplot(111, frameon=False)
	plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
	plt.grid(False)
	plt.xlabel("number of samples, $N$")
	plt.ylabel(r"$\tau$ estimates")
	for i in range(nrows):
		for j in range(ncols):
			dim = i*ncols + j
			taus = ax[i,j].loglog(N, tau[:,dim], "o-", label="new")
			line = ax[i,j].plot(N, N / 50.0, "--k", label=r"$\tau = N/50$")
	fname = "../results/"+objname+"/autocorr" + filename + ".png"
	fig.savefig(fname, format='png')
'''
