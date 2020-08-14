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
import pandas as pd
import emcee
import sys
import mm_likelihood


#chain = (nwalkers, nlink, ndim)
def plots(sampler, parameters, objname, fit_scale, float_names, obsdf, runprops, geo_obj_pos, mm_make_geo_pos):
			# Here parameters is whatever file/object will have the run params
	undo_ecc_aop = np.zeros(runprops.get('numobjects')-1)
	undo_ecc_aop[:] = False
	ecc_aop_index = np.zeros(runprops.get('numobjects')*2)
	for i in range(runprops.get('numobjects')-1):
		if 'ecc_'+str(i+2) in float_names and 'aop_'+str(i+2) in float_names:
			undo_ecc_aop[i] = True
			ecc_aop_index[2*i] = [float_names.index('ecc_'+str(i+2))]
			ecc_aop_index[2*i+1] = [float_names.index('aop_'+str(i+2))]
	flatchain = sampler.get_chain(flat = True)
	fit = []

	for i in fit_scale.columns:
		if i[0] in float_names:
			val = fit_scale.loc[0, i]
			fit.append(val)
            
	chain = sampler.get_chain(flat = False)            
	numparams = chain.shape[2]
	numwalkers = chain.shape[1]
	numgens = chain.shape[0]

	#First fit the flatchain with the fit parameters    
	fchain = np.zeros((numgens*numwalkers,numparams))    
	for i in range(numgens*numwalkers):
		row = np.zeros(numparams)        
		for j in range(numparams):
			val = flatchain[i][j]*fit[j]
			row[j] = val
		for i in range(runprops.get('numobjects')-1):           
			if undo_ecc_aop[i]:    
				aop_new = row[ecc_aop_index[i*2+1]]
				ecc_new = row[ecc_aop_index[i*2]]
				row[ecc_aop_index[i*2+1]] = np.arctan(aop_new/ecc_new)*180/np.pi
				row[ecc_aop_index[i*2]] = ecc_new/np.sin(aop_new)
		fchain[i] = row

	flatchain = np.array(fchain)


	#Now fit the chain 
	cchain = np.zeros((numgens,numwalkers, numparams))    
	for i in range(numgens):
		for j in range(numwalkers):
			row = []
			for k in range(numparams):
				val = chain[i][j][k]*fit[k]
				cchain[i][j][k] = val
			for i in range(runprops.get('numobjects')-1):
				if undo_ecc_aop[i]:    
					aop_new = cchain[i][j][ecc_aop_index[i*2+1]]
					ecc_new = cchain[i][j][ecc_aop_index[i*2]]
					cchain[i][j][ecc_aop_index[i*2+1]] = np.arctan(aop_new/ecc_new)*180/np.pi
					cchain[i][j][ecc_aop_index[i*2]] = ecc_new/np.sin(aop_new)

	cchain = np.array(cchain)

	oldchain = chain
	chain = cchain
	names = []
	for i in float_names:
		names.append(i)


	# Make corner plot
	#plt.rc('text', usetex=True)
	fig = corner.corner(flatchain, bins = 40, show_titles = True, 
			    plot_datapoints = False, color = "blue", fill_contours = True,
			    title_fmt = ".4e")
	fig.tight_layout(pad = 1.08, h_pad = 0, w_pad = 0)
	for ax in fig.get_axes():
		ax.tick_params(axis = "both", labelsize = 20, pad = 0.5)
	fname = "../runs/"+objname+"_"+runprops.get("date")+"/corner.pdf"       
	fig.savefig(fname, format = 'pdf')
	plt.close("all")
	#plt.rc('text', usetex=False)
	

	# Now make the walker plots
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
	sigsdf = pd.DataFrame(columns = ['-3sigma','-2sigma','-1sigma','median','1sigma','2sigma','3sigma', 'mean'], index = names)
	for i in range(len(flatchain[0])):        
		median = np.percentile(flatchain[:,i],50, axis = None)
		neg3sig= np.percentile(flatchain[:,i],0.37, axis = None)
		neg2sig = np.percentile(flatchain[:,i],2.275, axis = None)
		neg1sig = np.percentile(flatchain[:,i],15.866, axis = None)
		pos1sig = np.percentile(flatchain[:,i],84.134, axis = None)
		pos2sig = np.percentile(flatchain[:,i],97.724, axis = None)
		pos3sig = np.percentile(flatchain[:,i],99.63, axis = None)
		mean = np.mean(flatchain[:,i])
		sigsdf['-3sigma'].iloc[i] = neg3sig-median
		sigsdf['-2sigma'].iloc[i] = neg2sig-median
		sigsdf['-1sigma'].iloc[i] = neg1sig-median
		sigsdf['median'].iloc[i] = median
		sigsdf['1sigma'].iloc[i] = pos1sig-median
		sigsdf['2sigma'].iloc[i] = pos2sig-median
		sigsdf['3sigma'].iloc[i] = pos3sig-median
		sigsdf['mean'].iloc[i] = mean
	#if runprops.get('verbose'):
	print(sigsdf)
	filename = runprops.get('runs_folder') + '/sigsdf.csv'    
	sigsdf.to_csv(filename)
    
    
	for i in range(numparams):
		plt.figure(figsize = (9,9))
		plt.subplot(221)
		plt.hist(flatchain[:,i].flatten(), bins = 40, histtype = "step", color = "black")
		plt.subplot(223)
		plt.scatter(flatchain[:,i].flatten(), llhoods.flatten(),
			    c = np.mod(np.linspace(0,llhoods.size - 1, llhoods.size), numwalkers),
			    cmap = "nipy_spectral", edgecolors = "none")
		plt.xlabel(names[i])
		plt.ylabel("Log(L)")
		plt.subplot(224)
		plt.hist(llhoods.flatten(), bins = 40, orientation = "horizontal", 
			 histtype = "step", color = "black")
		plt.savefig("../runs/"+objname+"_"+runprops.get("date")+"/likelihood_" + names[i] + ".png")
		plt.close("all")

	# Residual plots
	llhoods = sampler.get_log_prob(flat = True)
	ind = np.argmin(llhoods)
	paramdf = flatchain[ind,:]

	chisquare_total, residuals = mm_likelihood.mm_chisquare(paramdf, obsdf, runprops, geo_obj_pos)

	colorcycle = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']
	objectnames = []
	for i in name_dict.values():
		objectnames.append(i)

	plt.figure()
	plt.Circle((0, 0), 1.0, color='black', fill=False)
	for i in range(1, nobjects):
		plt.scatter(residuals[2*(i-1)][:], residuals[2*(i-1)+1][:], c = colorcycle[i], label = objectnames[i], edgecolors = None)
	plt.xlabel("Delta Longitude")
	plt.ylabel("Delta Latitude")
	plt.axis("equal")
	plt.legend()
	plt.savefig("../runs/"+objname+"_"+runprops.get("date")+"/best_residuals.png")

	# Astrometry plots
	time_arr = obsdf['time'].values.flatten()
	tmin = time_arr.min()
	tmax = time_arr.max()
	fakeobsdf = obsdf.loc[[1,2],:]
	times = np.arange(tmin,tmax, 0.25)
	for i in range(len(times)):
		if i == 0 or i == 1:
			fakeobsdf.iloc[i,0] = times[i]
			# change row number?
		fakeobsdf = fakeobsdf.append(fakeobsdf.iloc[-1,:])
		fakeobsdf.iloc[-1,0] = times[i]
	geo_obj_pos = mm_make_geo_pos.mm_make_geo_pos("Haumea", times)

	llhoods = sampler.get_log_prob(flat = True)
	ind = np.argmin(llhoods)
	paramdf = flatchain[ind,:]

	Model_DeltaLong, Model_DeltaLat, fakeobsdf = mm_likelihood.mm_chisquare(paramdf, obsdf, runprops, geo_obj_pos, gensynth = True)

	modelx = np.empty((nobjects-1, obsdf.shape[0]))
	modely = np.empty((nobjects-1, obsdf.shape[0]))

	x = np.empty((nobjects-1, obsdf.shape[0]))
	xe = np.empty((nobjects-1, obsdf.shape[0]))
	y = np.empty((nobjects-1, obsdf.shape[0]))
	ye = np.empty((nobjects-1, obsdf.shape[0]))

	colorcycle = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']

	name_dict = runprops.get("names_dict")
	objectnames = []
	for i in name_dict.values():
		objectnames.append(i)

	fig = plt.figure()
	for i in range(1,nobjects):
		modelx[i-1,:] = Model_DeltaLong[i-1]
		modely[i-1,:] = Model_DeltaLat[i-1]

		x[i-1,:] = obsdf["DeltaLat_" + objectnames[i]].values
		xe[i-1,:] = obsdf["DeltaLat_" + objectnames[i] + "_err"].values
		y[i-1,:] = obsdf["DeltaLong_" + objectnames[i]].values
		ye[i-1,:] = obsdf["DeltaLong_" + objectnames[i] + "_err"].values

		plt.plot(modelx[i-1,:], modely[i-1,:], color = colorcycle[i], label = objectnames[i])
		plt.errorbar(x[i-1,:], y[i-1,:], xerr = xe[i-1,:], yerr = ye[i-1,:], fmt = "ko")

	plt.axis('equal')
	plt.xlabel("Delta Latitude")
	plt.ylabel("Delta Longitude")
	plt.legend()

	plt.savefig("../runs/"+objname+"_"+runprops.get("date")+"/best_astrometry.png")
	plt.close()



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
