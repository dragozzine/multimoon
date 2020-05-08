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
import numpy as np


def mm_analysis(sampler, parameters):
			# Here parameters is whatever file/object will have the run params
	flatchain = sampler.get_chain(flat = True)
	flatchain = np.transpose(flatchain)
	chain = sampler.get_chain(flat = False)
	# First start by converting the parameters into an array of strings
	# code here
	names = list(parameters)

	#fig = corner.corner(chain, bins = 40, labels = names, show_titles = True, 
	#		plot_datapoints = False, color = "blue", fill_contours = True,
	#		title_fmt = ".2f")
	#fig.tight_layout(pad = 1.08, h_pad = 0, w_pad = 0)
	#for ax in fig.get_axes():
	#	ax.tick_params(axis = "both", labelsize = 20, pad = 0.5)
	#fig.savefig(place to save the corner plot)

	
	# Now make the walker plots
	numsteps = chain.shape[0]
	numwalkers = chain.shape[1]
	numgens = chain.shape[2]
	print("numparams: ",numsteps)
	print("numwalkers: ",numwalkers)
	print("numgens: ", numgens)    
	print("names size: ",len(names))
	for i in range(len(names)):
		plt.figure()
		for j in range(numwalkers):
			plt.plot(np.reshape(chain[i,j,0:numgens],1))
		plt.ylabel(names[i])
		plt.xlabel("Generation")
		plt.savefig
		plt.close()

	# Likelihood plots??
	llhoods = sampler.get_log_prob(flat = True)
	for i in range(len(names)):
		plt.figure(figsize = (9,9))
		plt.subplot(221)
		plt.hist(chain[i,:,:].flatten(), bins = 40, histtype = "step", color = "black")
		plt.subplot(223)
		print(flatchain)
		print(len(flatchain[i,:].flatten()))
		print(len(llhoods.flatten()))
		plt.scatter(flatchain[i,:].flatten(), llhoods.flatten())
		plt.xlabel(names[i])
		plt.ylabel("Log(L)")
		plt.subplot(224)
		plt.hist(llhoods.flatten(), bins = 40, orientation = "horizontal", 
			 histtype = "step", color = "black")
		#plt.savefig()#place to save this)
		plt.close("all")