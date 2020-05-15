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



import corner
import matplotlib.pyplot as plt


def plots(sampler, parameters):
			# Here parameters is whatever file/object will have the run params
	flatchain = sampler.get_chain(flat = True)
	chain = sampler.get_chain(flat = False)
	# First start by converting the paramaters into an array of strings
	# code here
	names = list(paramaters)

	fig = corner.corner(chain, bins = 40, labels = names, show_titles = True, 
			    plot_datapoints = False, color = "blue", fill_contours = True,
			    title_fmt = ".4f")
	fig.tight_layout(pad = 1.08, h_pad = 0, w_pad = 0)
	for ax in fig.get_axes():
		ax.tick_params(axis = "both", labelsize = 20, pad = 0.5)
	#fig.savefig(place to save the corner plot)

	
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
		plt.savefig
		plt.close()

	# Likelihood plots??
	llhoods = sampler.get_log_prob(flat = True)
	for i in range(numparams):
		plt.figure(figsize = (9,9))
		plt.subplot(221)
		plt.hist(chain[i,:,:].flatten(), bins = 40, histtype = "step", color = "black")
		plt.subplot(223)
		plt.scatter(flatchain[i,:].flatten(), llhoods.flatten())
		plt.xlabel(names[i])
		plt.ylabel("Log(L)")
		plt.subplot(224)
		plt.hist(llhoods.flatten(), bins = 40, orientation = "horizontal", 
			 histtype = "step", color = "black")
		plt.savefig(#place to save this)
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

def autocorrelation(sampler, parameters):
	# Getting chain for the first parameter to calculate different values
	chain = sampler.get_chain()
	
	nwalkers = sampler.nwalkers
	ndims = sampler.ndim
	nsteps = chain.shape[0]
	
	# Converting parameters df to array of names
	names = list(parameters)

	# Calculating values to calculate tau for
	# This chould be changed eventually
	N = np.exp(np.linspace(np.log(100), np.log(nsteps), 10)).astype(int)

	# Setting up array for tau estimates
	tau = np.empty( (len(N), ndims) )

	# Calculating tau for each value in N for each parameter
	for i in range(ndims):
		thischain = chain[:,:,i].T
		for j, n in enumerate(N):
			tau[j,i] = autocorr_new(chain[:, :n])

	# Setting up to plot curves for tau in a grid
	x = 3
	y = ndims
	nrows = 0
	ncols = 3
	while x <= y:
		y = y - x
		nrows += 1

	# Plotting
	fig, ax = plt.subplots(nrows = nrows, ncols = ncols, sharey=True, 
			       gridspec_kw={'wspace': 0},
			       figsize = (6.4*(ncols-1),4.8*(nrows)), 
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
			taus = ax[i,j].loglog(N, new[:,dim], "o-", label="new")
			line = ax[i,j].plot(N, N / 50.0, "--k", label=r"$\tau = N/50$")
			ax[i,j].text(N[0], 100, names[dim])
			if i == 0 and j == 0:
				 plt.legend([l3],[r"$\tau = N/50$"],fontsize = 14,
					    loc="lower right")

	#fig.savefig("autocorr.png")
	# Outputting the emcee calculated autocorrelation time as an additional check
	print(
		"Mean autocorrelation time: {0:.3f} steps".format(
			np.mean(sampler.get_autocorr_time())
		)
	)
	print("Estimated autocorrelation time for each parameter:")
	print(sampler.get_autocorr_time())
