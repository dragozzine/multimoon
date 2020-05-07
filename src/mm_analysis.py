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


import corner
import matplotlib.pyplot as plt


def mm_analysis(sampler, parameters):
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
	numparams = chain.shape[0]
	numwalkers = chain.shape[1]
	numgens = chain.shape[2]
	for i in range(numparams):
		plt.figure()
		for j in range(numwalkers):
			plt.plot(np.reshape(chain[i,j,0:numgens))
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


