import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import corner
import numpy as np
import pandas as pd
import emcee
import sys
import os
import mm_likelihood
from astropy.time import Time
import commentjson as json
import mm_param
import mm_make_geo_pos
from mm_SPINNY.spinny_plots import spinny_plot
from mm_SPINNY.spinny_generate import *
from mm_SPINNY.spinny_vector import *
from mm_SPINNY.spinny_nosun import *
from mm_SPINNY.mm_vpython import *
from mm_SPINNY.keplerian import *

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data
def save(s_df,names):

    save_yn = "N"
    if save_yn=="Y":
        print("Generating .csv...")
        t_current = ctime().replace(" ","_")
        filename = names[1]+"_SPINNY_"+t_current+".csv"
        s_df.to_csv("output/"+filename)
        print("SPINNY data saved to the output file as "+filename)
        plot_q(s_df, names)
    elif save_yn == "N":
        print("")
        plot_q(s_df, names)
    else:
        print("")
        print('Invalid Response.')
        return save(s_df,names)   
    
    
#chain = (nwalkers, nlink, ndim)

def plots(sampler, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names):
	# Here total_df_names is whatever file/object will have the run params
	objname = runprops.get("objectname")
	nthinning = runprops.get("nthinning")
	undo_ecc_aop = np.zeros(runprops.get('numobjects')-1)
	undo_ecc_aop[:] = False
	ecc_aop_index = np.zeros((runprops.get('numobjects')-1)*2)
	undo_inc_lan = np.zeros(runprops.get('numobjects')-1)
	undo_inc_lan[:] = False
	undo_spin = np.zeros(runprops.get('numobjects'))
	undo_spin[:] = False
	spin_index = np.zeros((runprops.get('numobjects')-1)*2)
	inc_lan_index = np.zeros((runprops.get('numobjects')-1)*2)
	undo_lambda = np.zeros(runprops.get('numobjects')-1)
	undo_lambda[:] = False
	lambda_index = np.zeros((runprops.get('numobjects')-1)*2)
	undo_pomega = np.zeros(runprops.get('numobjects')-1)
	undo_pomega[:] = False
	pomega_index = np.zeros((runprops.get('numobjects')-1)*2)
	undo_masses = np.zeros(2)    
	undo_masses[:] = False
	masses_index = np.zeros(runprops.get('numobjects'))
	#print(os.getcwd())    
    
	for i in range(runprops.get('numobjects')-1):
		if 'ecc_'+str(i+2) in float_names and 'aop_'+str(i+2) in float_names:
			undo_ecc_aop[i] = True
			ecc_aop_index[2*i] = float_names.index('ecc_'+str(i+2))
			ecc_aop_index[2*i+1] = float_names.index('aop_'+str(i+2))
		if 'inc_'+str(i+2) in float_names and 'lan_'+str(i+2) in float_names:
			undo_inc_lan[i] = True
			inc_lan_index[2*i] = float_names.index('inc_'+str(i+2))
			inc_lan_index[2*i+1] = float_names.index('lan_'+str(i+2))
		if 'spinc_'+str(i+2) in float_names and 'splan_'+str(i+2) in float_names:
			undo_spin[i] = True
			spin_index[2*i] = float_names.index('spinc_'+str(i+2))
			spin_index[2*i+1] = float_names.index('splan_'+str(i+2))
		if 'mea_'+str(i+2) in float_names and 'aop_'+str(i+2) in float_names:
			undo_lambda[i] = True
			lambda_index[2*i] = float_names.index('mea_'+str(i+2))
			lambda_index[2*i+1] = float_names.index('aop_'+str(i+2))
		if 'aop_'+str(i+2) in float_names and 'lan_'+str(i+2) in float_names:
			undo_pomega[i] = True
			pomega_index[2*i] = float_names.index('aop_'+str(i+2))
			pomega_index[2*i+1] = float_names.index('lan_'+str(i+2))
	if 'mass_1' in float_names and 'mass_2' in float_names:
		if 'mass_3' in float_names and runprops.get('numobjects') > 2:        
			undo_masses[1] = True
			masses_index[0] = float_names.index('mass_1')
			masses_index[1] = float_names.index('mass_2')
			masses_index[2] = float_names.index('mass_3')
		else:        
			undo_masses[0] = True
			masses_index[0] = float_names.index('mass_1')
			masses_index[1] = float_names.index('mass_2')

	burnin = int(runprops.get('nburnin'))
	clusterburn = int(runprops.get('clustering_burnin'))
	thin_plots = runprops.get('thin_plots')    
	chain = sampler.get_chain(discard=int(burnin+clusterburn),flat = False, thin=thin_plots)
	fit = []

	for i in fit_scale.columns:
		name = i
		if type(name) != str:
			name = name[0]
		if name in float_names:
			val = fit_scale.loc[0, i]
			fit.append(val)

	# Getting final values for the shape of the chain
	shortchain = sampler.get_chain(discard=int(burnin+clusterburn),flat = False, thin=thin_plots)
	numparams = shortchain.shape[2]
	numwalkers = shortchain.shape[1]
	numgens = shortchain.shape[0]
	del shortchain	
 
	# Take chain "fit" values and make them into real values
	for i in range(numparams):
		chain[:,:,i] = chain[:,:,i]*fit[i]

	fitparam_chain = np.zeros((1,numwalkers,numgens))
	print(fitparam_chain.shape)    
	fitparam_names = []    
	# Now de-transform the chain
	print("Starting un transformations")
	if runprops.get("transform"):
		for b in range(runprops.get('numobjects')-1):
			if undo_ecc_aop[b]:
				aop_new = chain[:,:,int(ecc_aop_index[b*2+1])]
				ecc_new = chain[:,:,int(ecc_aop_index[b*2])]
				print(aop_new.T.shape, np.array([aop_new.T]).shape)
				fitparam_chain = np.concatenate((fitparam_chain, np.array([aop_new.T])),axis=0)
				fitparam_chain = np.concatenate((fitparam_chain, np.array([ecc_new.T])),axis=0)                
				fitparam_names.append('aop_new')
				fitparam_names.append('ecc_new')
				pomega = (np.arctan2(ecc_new,aop_new)*180/np.pi)%360
				chain[:,:,int(ecc_aop_index[b*2+1])] = pomega
				chain[:,:,int(ecc_aop_index[b*2])] = ecc_new/np.sin(pomega/180*np.pi)
			if undo_inc_lan[b]:
				inc_new = chain[:,:,int(inc_lan_index[b*2])]
				lan_new = chain[:,:,int(inc_lan_index[b*2+1])]
				fitparam_chain = np.concatenate((fitparam_chain, np.array([inc_new.T])),axis=0)
				fitparam_chain = np.concatenate((fitparam_chain, np.array([lan_new.T])),axis=0)
				fitparam_names.append('inc_new')
				fitparam_names.append('lan_new')
				lan = (np.arctan2(inc_new,lan_new)*180/np.pi)%360
				chain[:,:,int(inc_lan_index[b*2+1])] = lan
				inc = (np.arctan2(inc_new,np.sin(lan*np.pi/180))*2*180/np.pi)%180
				chain[:,:,int(inc_lan_index[b*2])] = inc
			if undo_lambda[b]:
				mea_new = chain[:,:,int(lambda_index[b*2])]
				pomega = chain[:,:,int(lambda_index[b*2+1])]
				fitparam_chain = np.concatenate((fitparam_chain, np.array([mea_new.T])),axis=0)
				fitparam_chain = np.concatenate((fitparam_chain, np.array([pomega.T])),axis=0)
				fitparam_names.append('mea_new')
				fitparam_names.append('pomega')
				mea = (mea_new-pomega)%360
				chain[:,:,int(lambda_index[b*2])] = mea
			if undo_pomega[b]:
				lan = chain[:,:,int(pomega_index[b*2+1])]
				pomega = chain[:,:,int(pomega_index[b*2])]
				aop = (pomega-lan)%360
				chain[:,:,int(pomega_index[b*2])] = aop
		for b in range(runprops.get('numobjects')):
			if undo_spin[b]:
				spinc_new = chain[:,:,int(spin_index[b*2])]
				splan_new = chain[:,:,int(spin_index[b*2+1])]
				fitparam_chain = np.concatenate((fitparam_chain, np.array([spinc_new.T])),axis=0)
				fitparam_chain = np.concatenate((fitparam_chain, np.array([splan_new.T])),axis=0)
				fitparam_names.append('sp_p')
				fitparam_names.append('sp_q')
				splan = (np.arctan2(spinc_new,splan_new)*180/np.pi)%360
				chain[:,:,int(spin_index[b*2+1])] = lan
				spinc = (np.arctan2(spinc_new,np.sin(splan*np.pi/180))*2*180/np.pi)%180
				chain[:,:,int(spin_index[b*2])] = spinc
		if undo_masses[0]:
			mass_1 = chain[:,:,int(masses_index[0])]
			mass_2 = chain[:,:,int(masses_index[1])]
			fitparam_chain = np.concatenate((fitparam_chain, np.array([mass_2.T])),axis=0)
			fitparam_names.append('mass1+2')
			chain[:,:,int(masses_index[1])] = mass_2-mass_1
		elif undo_masses[1]:
			mass_1 = chain[:,:,int(masses_index[0])]
			mass_2 = chain[:,:,int(masses_index[1])]
			mass_3 = chain[:,:,int(masses_index[2])]
			fitparam_chain = np.concatenate((fitparam_chain, np.array([mass_2.T])),axis=0)
			fitparam_chain = np.concatenate((fitparam_chain, np.array([mass_3.T])),axis=0)
			fitparam_names.append('mass1+2')
			fitparam_names.append('mass1+2+3')
			chain[:,:,int(masses_index[2])] = (mass_3-mass_2)/(10**18) 
			chain[:,:,int(masses_index[1])] = (mass_2-mass_1)/(10**18)
			chain[:,:,int(masses_index[0])] = (mass_1)/(10**18)

	fitparam_chain = np.delete(fitparam_chain,0,0)
	fitparam_chain = fitparam_chain.T    
        
	print("Un transforming done")

	# Cutting up chain
	full_chain = np.copy(chain)
	#chain = chain[int(burnin+clusterburn + thin_plots - 1) :: thin_plots]
	print(chain.shape)

	# Flattening the chain based on method in emcee
	s = list(chain.shape[1:])
	s[0] = np.prod(chain.shape[:2])
	s2 = list(fitparam_chain.shape[1:])
	s2[0] = np.prod(fitparam_chain.shape[:2])
	print(s2, fitparam_chain.shape)    
	flatchain = chain.reshape(s)
	fitparam_chain = fitparam_chain.reshape(s2)    
	print(flatchain.shape, fitparam_chain.shape)

	# Getting parameter names
	names = []
	#print(float_names)    
	for i in float_names:
		names.append(i)

	# Getting log likelihood posterior values for use throughout
	llhoods = sampler.get_log_prob(discard=int(burnin+clusterburn),flat = True, thin=thin_plots)
	print(llhoods.shape)    
	ind = np.argmax(llhoods)
	params = flatchain[ind,:].flatten()

	# Making latex labels for values
	latexnames = names.copy()
	for i in range(len(latexnames)):
		if "mass" in latexnames[i]:
			latexnames[i] = "$"+latexnames[i].replace("mass","m")+"$ ($10^{18}$ kg)"
			latexnames[i] = fr"{latexnames[i]}"
		if "sma" in latexnames[i]:
			latexnames[i] = "$"+latexnames[i].replace("sma","a")+"$ (km)"
			latexnames[i] = fr"{latexnames[i]}"
		if "ecc" in latexnames[i]:
			latexnames[i] = "$"+latexnames[i].replace("ecc","e")+"$"
			latexnames[i] = fr"{latexnames[i]}"
		if "aop" in latexnames[i]:
			latexnames[i] = "$"+latexnames[i].replace("aop","\omega")+"$ ($^{\circ}$)"
			latexnames[i] = fr"{latexnames[i]}"
		if "inc" in latexnames[i]:
			latexnames[i] = "$"+latexnames[i].replace("inc","i")+"$ ($^{\circ}$)"
			latexnames[i] = fr"{latexnames[i]}"
		if "lan" in latexnames[i]:
			latexnames[i] = "$"+latexnames[i].replace("lan","\Omega")+"$ ($^{\circ}$)"
			latexnames[i] = fr"{latexnames[i]}"
		if "mea" in latexnames[i]:
			latexnames[i] = "$"+latexnames[i].replace("mea","M")+"$ ($^{\circ}$)"
			latexnames[i] = fr"{latexnames[i]}"
		if "j2r2" in latexnames[i]:
			latexnames[i] = "$"+latexnames[i].replace("j2r2","J_2R^2")+"$ (km$^2$)"
			latexnames[i] = fr"{latexnames[i]}"
		if "spinc" in latexnames[i]:
			latexnames[i] = "$"+latexnames[i].replace("spinc","i^{spin}")+"$ ($^{\circ}$)"
			latexnames[i] = fr"{latexnames[i]}"
		if "splan" in latexnames[i]:
			latexnames[i] = "$"+latexnames[i].replace("splan","\Omega^{spin}")+"$ ($^{\circ}$)"
			latexnames[i] = fr"{latexnames[i]}"
		if "c22r2" in latexnames[i]:
			latexnames[i] = "$"+latexnames[i].replace("c22r2","C_{22}R^2")+"$ (km$^2$)"
			latexnames[i] = fr"{latexnames[i]}"
		if "spaop" in latexnames[i]:
			latexnames[i] = "$"+latexnames[i].replace("spaop","\omega^{spin}")+"$ ($^{\circ}$)"
			latexnames[i] = fr"{latexnames[i]}"
		if "sprate" in latexnames[i]:
			latexnames[i] = "$"+latexnames[i].replace("sprate","")+"$ (rad s$^{-1}$)"
			latexnames[i] = fr"{latexnames[i]}"

	# Make corner plot
	#plt.rc('text', usetex=True)
	fig = 0
	fig = corner.corner(flatchain, labels = latexnames, bins = 40, show_titles = True, 
			    plot_datapoints = False, color = "blue", fill_contours = True,
			    title_fmt = ".3f", truths = params, label_kwargs=dict(fontsize=20))
	fig.tight_layout(pad = 1.08, h_pad = -0.4, w_pad = -0.4)
	for ax in fig.get_axes():
		ax.tick_params(axis = "both", labelsize = 12, pad = 0.0)
	fname = "corner.pdf"       
	fig.savefig(fname, format = 'pdf')
	plt.close("all")

	# Making corner plots with derived parameters
	dnames = names.copy()
	dfchain = flatchain.copy()
	dtif = runprops.get("dynamicstoincludeflags")

	# Orbital periods for each satellite
	for i in range(1,runprops.get('numobjects')):
		print(names)        
		a_index = [n for n, l in enumerate(names) if l.startswith('sma_')][0]
		m_index = [n for n, l in enumerate(names) if l.startswith('mass_')][0]
		mp_index = [n for n, l in enumerate(names) if l.startswith('mass_1')][0]

		a_arr = flatchain[:,a_index]
		m_arr = flatchain[:,m_index]
		mp_arr = flatchain[:,mp_index]

		period = 2*np.pi*np.sqrt(a_arr**3/(6.674e-20*(m_arr + mp_arr)*10**18))/3600.0/24.0

		dnames = np.append(dnames, ["period_" + str(i+1)])
		dfchain = np.concatenate((dfchain, np.array([period]).T), axis = 1)

	# Spin period
	fixedparams = runprops.get("float_dict")
	if dtif[0] == "1" or dtif[0] == "2":
		if fixedparams["sprate_1"] == 1:
			sprate1_index = [n for n, l in enumerate(names) if l.startswith('sprate_1')][0]
			sprate1_arr = flatchain[:,sprate1_index]
			spinperiod = (2*np.pi/sprate1_arr)/3600.0

			dnames = np.append(dnames, ["spin_period_1"])
			dfchain = np.concatenate((dfchain, np.array([spinperiod]).T), axis = 1)

	# Satellite-spin mutual inclination
	if dtif[0] == "1" or dtif[0] == "2":
		for i in range(1,runprops.get('numobjects')):
			print(names)            
			spinc1_index = [n for n, l in enumerate(names) if l.startswith('spinc_1')][0]
			splan1_index = [n for n, l in enumerate(names) if l.startswith('splan_1')][0]
			inc_index = [n for n, l in enumerate(names) if l.startswith('inc_'+str(i+1))][0]
			lan_index = [n for n, l in enumerate(names) if l.startswith('lan_'+str(i+1))][0]

			spinc1_arr = np.deg2rad(flatchain[:,spinc1_index])
			splan1_arr = np.deg2rad(flatchain[:,splan1_index])
			inc_arr = np.deg2rad(flatchain[:,inc_index])
			lan_arr = np.deg2rad(flatchain[:,lan_index])

			mutualinc = np.arccos( np.cos(spinc1_arr)*np.cos(inc_arr) + np.sin(spinc1_arr)*np.sin(inc_arr)*np.cos(splan1_arr - lan_arr) )
			mutualinc = np.rad2deg(mutualinc)

			dnames = np.append(dnames, ["sat-spin inc_" + str(i+1)])
			dfchain = np.concatenate((dfchain, np.array([mutualinc]).T), axis = 1)

	# Satellite-Satellite mutual inclination (this only works right now for 2 moons/satellites)
	if runprops.get("numobjects") > 2:
		inc2_index = [n for n, l in enumerate(names) if l.startswith('inc_2')][0]
		inc3_index = [n for n, l in enumerate(names) if l.startswith('inc_3')][0]
		lan2_index = [n for n, l in enumerate(names) if l.startswith('lan_2')][0]
		lan3_index = [n for n, l in enumerate(names) if l.startswith('lan_3')][0]

		inc2_arr = np.deg2rad(flatchain[:,inc2_index])
		lan2_arr = np.deg2rad(flatchain[:,lan2_index])
		inc3_arr = np.deg2rad(flatchain[:,inc3_index])
		lan3_arr = np.deg2rad(flatchain[:,lan3_index])

		mutualinc = np.arccos( np.cos(inc2_arr)*np.cos(inc3_arr) + np.sin(inc2_arr)*np.sin(inc3_arr)*np.cos(lan2_arr - lan3_arr) )
		mutualinc = np.rad2deg(mutualinc)

		dnames = np.append(dnames, ["sat-sat inc"])
		dfchain = np.concatenate((dfchain, np.array([mutualinc]).T), axis = 1)

# Creating corner+derived plot
	fig = corner.corner(dfchain, labels = dnames, bins = 40, show_titles = True, 
			    plot_datapoints = False, color = "blue", fill_contours = True,
			    title_fmt = ".4f", truths = dfchain[ind,:].flatten())
	#fig.tight_layout(pad = 1.08, h_pad = 0, w_pad = 0)
	#for ax in fig.get_axes():
	#	ax.tick_params(axis = "both", labelsize = 20, pad = 0.5)
	fname = "corner+derived.pdf"       
	fig.savefig(fname, format = 'pdf')
	plt.close("all")


	# Creating corner_fitparams plot
	fitflatchain = sampler.get_chain(discard=int(burnin+clusterburn),flat = True, thin=thin_plots)

	fig = corner.corner(fitflatchain, bins = 40, show_titles = True, 
			    plot_datapoints = False, color = "blue", fill_contours = True,
			    title_fmt = ".6f", truths = fitflatchain[ind,:].flatten())
	#fig.tight_layout(pad = 1.08, h_pad = 0, w_pad = 0)
	#for ax in fig.get_axes():
	#	ax.tick_params(axis = "both", labelsize = 20, pad = 0.5)
	fname = "corner_fitparams.pdf"       
	fig.savefig(fname, format = 'pdf')
	plt.close("all")
	
	del fitflatchain


	#plt.rc('text', usetex=False)
	

	# Now make the walker plots
	from matplotlib.backends.backend_pdf import PdfPages

	walkerpdf = PdfPages("walkers.pdf")

	for i in range(numparams):
		plt.figure(dpi = 50)
		for j in range(numwalkers):
			plt.plot(np.reshape(chain[0:numgens,j,i], numgens), alpha=0.2)
		plt.ylabel(names[i])
		plt.xlabel("Generation")
		#plt.savefig(runprops.get('results_folder')+"/walker_"+names[i]+".png")
		walkerpdf.attach_note(names[i])
		walkerpdf.savefig()
		#plt.close()

	walkerpdf.close()
	plt.close("all")
    
	fullwalkerpdf = PdfPages("walkers_full.pdf")
	backend = emcee.backends.HDFBackend('chain.h5')    

#	full_chain = sampler.get_chain(discard=0, flat = False)  
	fullgens = full_chain.shape[0]
	#print(fullgens)
	for i in range(numparams):
		plt.figure(dpi = 50)
		for j in range(numwalkers):
			plt.plot(np.reshape(full_chain[0:fullgens,j,i], fullgens), alpha=0.2)
		plt.axvline(x=burnin)
		plt.axvline(x=(clusterburn+burnin))
		plt.ylabel(names[i])
		plt.xlabel("Generation")
		#plt.savefig(runprops.get('results_folder')+"/walker_"+names[i]+".png")
		fullwalkerpdf.attach_note(names[i])
		fullwalkerpdf.savefig()
		#plt.close()

	fullwalkerpdf.close()
	plt.close("all")

	# Figuring out the distributions of total_df_names
	#old_fchain = sampler.get_chain(flat=True)
	llhoods = sampler.get_log_prob(discard=int(burnin+clusterburn),flat = True, thin=thin_plots)
	sigsdf = pd.DataFrame(columns = ['-3sigma','-2sigma','-1sigma','median','1sigma','2sigma','3sigma', 'mean', 'best fit'], index = dnames)
	j = 0
	for i in range(len(dfchain[0])):
		num = dfchain[:,i]
#		if i>len(names):
# 			
		median = np.percentile(num,50, axis = None)
		neg3sig= np.percentile(num,0.37, axis = None)
		neg2sig = np.percentile(num,2.275, axis = None)
		neg1sig = np.percentile(num,15.866, axis = None)
		pos1sig = np.percentile(num,84.134, axis = None)
		pos2sig = np.percentile(num,97.724, axis = None)
		pos3sig = np.percentile(num,99.63, axis = None)
		mean = np.mean(num)
		bestfit = dfchain[ind,:].flatten()[i]
		sigsdf['-3sigma'].iloc[i] = neg3sig-median
		sigsdf['-2sigma'].iloc[i] = neg2sig-median
		sigsdf['-1sigma'].iloc[i] = neg1sig-median
		sigsdf['median'].iloc[i] = median
		sigsdf['1sigma'].iloc[i] = pos1sig-median
		sigsdf['2sigma'].iloc[i] = pos2sig-median
		sigsdf['3sigma'].iloc[i] = pos3sig-median
		sigsdf['mean'].iloc[i] = mean
		sigsdf['best fit'].iloc[i] = bestfit
	#if runprops.get('verbose'):
	print(sigsdf)
	filename = 'sigsdf.csv'    
	sigsdf.to_csv(filename)
    
	# Likelihood plots    
	likelihoodspdf = PdfPages("likelihoods.pdf")
	ylimmin = np.percentile(llhoods.flatten(), 1)
	ylimmax = llhoods.flatten().max() + 1
	print(chain.shape,flatchain.shape, llhoods.shape)
	dfparams = dfchain.shape[1]
	for i in range(numparams):
		plt.figure(figsize = (9,9))
		plt.subplot(221)
		plt.hist(flatchain[:,i].flatten(), bins = 40, histtype = "step", color = "black")
		plt.subplot(223)
		plt.scatter(flatchain[:,i].flatten(), llhoods.flatten(),
			    c = np.mod(np.linspace(0,llhoods.size - 1, llhoods.size), numwalkers),
			    cmap = "nipy_spectral", edgecolors = "none", rasterized=True, alpha=0.1)
		plt.xlabel(dnames[i])
		plt.ylabel("Log(L)")
		plt.ylim(ylimmin, ylimmax)
		plt.subplot(224)
		llflat = llhoods.flatten()
		plt.hist(llflat[np.isfinite(llflat)], bins = 40, orientation = "horizontal", 
			 histtype = "step", color = "black")
		plt.ylim(ylimmin, ylimmax)
		likelihoodspdf.attach_note(dnames[i])
		likelihoodspdf.savefig()
		#plt.savefig(runprops.get("results_folder")+"/likelihood_" + names[i] + ".png")
        
	for i in range(numparams,dfchain.shape[1]):
		plt.figure(figsize = (9,9))
		plt.subplot(221)
		plt.hist(dfchain[:,i].flatten(), bins = 40, histtype = "step", color = "black")
		plt.subplot(223)
		plt.scatter(dfchain[:,i].flatten(), llhoods.flatten(),
			    c = np.mod(np.linspace(0,llhoods.size - 1, llhoods.size), numwalkers),
			    cmap = "nipy_spectral", edgecolors = "none", rasterized=True, alpha=0.1)
		plt.xlabel(dnames[i])
		plt.ylabel("Log(L)")
		plt.ylim(ylimmin, ylimmax)
		plt.subplot(224)
		llflat = llhoods.flatten()
		plt.hist(llflat[np.isfinite(llflat)], bins = 40, orientation = "horizontal", 
			 histtype = "step", color = "black")
		plt.ylim(ylimmin, ylimmax)
		likelihoodspdf.attach_note(dnames[i])
		likelihoodspdf.savefig()
		#plt.savefig(runprops.get("results_folder")+"/likelihood_" + names[i] + ".png")
          
	for i in range(fitparam_chain.shape[1]):
		plt.figure(figsize = (9,9))
		plt.subplot(221)
		plt.hist(fitparam_chain[:,i].flatten(), bins = 40, histtype = "step", color = "black")
		plt.subplot(223)
		plt.scatter(fitparam_chain[:,i].flatten(), llhoods.flatten(),
			    c = np.mod(np.linspace(0,llhoods.size - 1, llhoods.size), numwalkers),
			    cmap = "nipy_spectral", edgecolors = "none", rasterized=True, alpha=0.1)
		plt.xlabel(fitparam_names[i])
		plt.ylabel("Log(L)")
		plt.ylim(ylimmin, ylimmax)
		plt.subplot(224)
		llflat = llhoods.flatten()
		plt.hist(llflat[np.isfinite(llflat)], bins = 40, orientation = "horizontal", 
			 histtype = "step", color = "black")
		plt.ylim(ylimmin, ylimmax)
		likelihoodspdf.attach_note(fitparam_names[i])
		likelihoodspdf.savefig()
		#plt.savefig(runprops.get("results_folder")+"/likelihood_" + names[i] + ".png")

	likelihoodspdf.close()
	plt.close("all")

	# Residual plots
	flatchain = sampler.get_chain(flat = True)
	nobjects = runprops.get('numobjects')
	llhoods = sampler.get_log_prob(flat = True)
	ind = np.argmax(llhoods)
	params = flatchain[ind,:].flatten()
    
	paraminput = []

	for i in params:
		paraminput.append(i)
        

	paramnames = names
	name_dict = runprops.get("names_dict")

	objectnames = []
	for i in range(runprops.get('numobjects')):
		objectnames.append(name_dict.get('name_'+str(i)))
    
	for values,keys in name_dict.items():
		for j in range(runprops.get('numobjects')):
			if str(j+1) in keys or 'offset' in keys: 
				paraminput.append(values)
				paramnames.append(keys)

	#print(paraminput)
	#print(paramnames)
	names_dict = runprops.get("names_dict")    
	paramdf,fit_params = mm_param.from_fit_array_to_param_df(paraminput, paramnames, fixed_df, total_df_names, fit_scale, names_dict, runprops)

	#paramdf = pd.DataFrame(paraminput).transpose()
	#paramdf.columns = paramnames


	#print(paramdf)
#Currently this function call sends an error in the case of leaving any necessary value floating, since paramdf will be incomplete 
	#chisquare_total, residuals = mm_likelihood.mm_chisquare(paramdf, obsdf, runprops, geo_obj_pos)
	best_likelihoods = pd.read_csv('best_likelihoods.csv')
	residuals = []
	#print(best_likelihoods, best_likelihoods.iloc[-(i+1)])    
	for i in range(runprops.get('numobjects')):    
		residuals.insert(0,best_likelihoods.iloc[-(i+1)][-1])
		residuals.insert(0,best_likelihoods.iloc[-(i+2)][-1])
	#print(residuals)        
	#print(chisquare_total, residuals)

	colorcycle = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']

	xvals1 = np.linspace(-1.0,1.0,num=1000)
	xvals2 = np.linspace(-2.0,2.0,num=1000)
	xvals3 = np.linspace(-3.0,3.0,num=1000)
	circle1 = np.sqrt(1 - xvals1**2)
	circle2 = np.sqrt(4 - xvals2**2)
	circle3 = np.sqrt(9 - xvals3**2)

	plt.figure()
	plt.plot(xvals1, circle1, color = "black")
	plt.plot(xvals1,-circle1, color = "black")
	plt.plot(xvals2, circle2, color = "black", alpha = 0.5)
	plt.plot(xvals2,-circle2, color = "black", alpha = 0.5)
	plt.plot(xvals3, circle3, color = "black", alpha = 0.25)
	plt.plot(xvals3,-circle3, color = "black", alpha = 0.25)
	print(nobjects, np.array(residuals).shape)    
	for i in range(1, nobjects):
		plt.scatter(residuals[2*(i-1)][:], residuals[2*(i-1)+1][:], c = colorcycle[i], label = objectnames[i], edgecolors = None)
	plt.xlabel("Delta Longitude")
	plt.ylabel("Delta Latitude")
	plt.axis("equal")
	plt.legend()
	plt.savefig("best_residuals.pdf", format = "pdf")

	# Astrometry plots
	time_arr = obsdf['time'].values.flatten()
	tmin = time_arr.min()
	tmax = time_arr.max()

	converttimes = [tmin,tmax]
	t = Time(converttimes, format = 'jd')

	timesdic = {'start': t.isot[0], 'stop': t.isot[1], 'step': '6h'}
    
	#geo_obj_pos = mm_make_geo_pos.mm_make_geo_pos(objname, timesdic, runprops, True)
	geo_obj_pos = pd.read_csv('geocentric_'+objname+'_position_analysis.csv')

	times = geo_obj_pos.values[:,0].flatten()

	fakeobsdf = obsdf.loc[[0,1],:]
	for i in range(len(times)):
		if i == 0 or i == 1:
			fakeobsdf.iloc[i,0] = times[i]
			# change row number?
		fakeobsdf = fakeobsdf.append(fakeobsdf.iloc[-1,:])
		fakeobsdf['time'].iloc[-1] = times[i]
	fakeobsdf = fakeobsdf.iloc[2:]

	names = list(runprops.get('names_dict').values())
	vals = ['mass','axis','j2r2','c22r2','sp_rate','sp_obl','sp_prc','sp_lon','sma','ecc','inc','lan','aop','mea']    
	sys_df = pd.DataFrame(columns=names,index=vals)
	sys_df[names[0]]['sma'] = 0
	sys_df[names[0]]['ecc'] = 0
	sys_df[names[0]]['inc'] = 0
	sys_df[names[0]]['lan'] = 0
	sys_df[names[0]]['aop'] = 0
	sys_df[names[0]]['mea'] = 0

	if runprops.get('includesun'):
		sys_df.loc['mass','Sun'] = paramdf['mass_0'][0]
		axis = list(runprops.get('axes_size').values())[0]
		j2 = 0
		c22 = 0
		sp_rate = 0
		sp_obl = 0
		sp_prc = 0
		sp_lon = 0      
		sys_df.loc['axis','Sun'] = axis      
		sys_df.loc['j2r2','Sun'] = j2     
		sys_df.loc['c22r2','Sun'] = c22    
		sys_df.loc['sp_rate','Sun'] = sp_rate
		sys_df.loc['sp_obl','Sun'] = sp_obl
		sys_df.loc['sp_prc','Sun'] = sp_prc
		sys_df.loc['sp_lon','Sun'] = sp_lon          
		sma = paramdf['sma_0'][0]
		ecc = paramdf['ecc_0'][0]
		inc = paramdf['inc_0'][0]
		lan = paramdf['lan_0'][0]
		aop = paramdf['aop_0'][0]
		mea = paramdf['mea_0'][0]
		sys_df.loc['sma','Sun'] = sma
		sys_df.loc['ecc','Sun'] = ecc
		sys_df.loc['inc','Sun'] = inc
		sys_df.loc['lan','Sun'] = lan
		sys_df.loc['aop','Sun'] = aop
		sys_df.loc['mea','Sun'] = mea
    
	for i in range(runprops.get('numobjects')):
		sys_df.loc['mass',names[i]] = paramdf['mass_'+str(i+1)][0]
		axis = list(runprops.get('axes_size').values())[i]
		if int(runprops.get('dynamicstoincludeflags')[i]) == 0:        
			sys_df.loc['axis',names[i]] = 0        
			sys_df.loc['j2r2',names[i]] = 0        
			sys_df.loc['c22r2',names[i]] = 0        
			sys_df.loc['sp_rate',names[i]] = 0        
			sys_df.loc['sp_obl',names[i]] = 0        
			sys_df.loc['sp_prc',names[i]] = 0        
			sys_df.loc['sp_lon',names[i]] = 0
		elif int(runprops.get('dynamicstoincludeflags')[i]) == 1:   
			j2 = paramdf['j2r2_'+str(i+1)][0]
			c22 = paramdf['c22r2_'+str(i+1)][0]
			sp_rate = paramdf['sprate_'+str(i+1)][0]
			sp_obl = paramdf['spinc_'+str(i+1)][0]
			sp_prc = paramdf['splan_'+str(i+1)][0]
			sp_lon = paramdf['spaop_'+str(i+1)][0]
			sys_df.loc['axis',names[i]] = axis
			sys_df.loc['j2r2',names[i]] = j2        
			sys_df.loc['c22r2',names[i]] = 0
			sys_df.loc['sp_rate',names[i]] = sp_rate
			sys_df.loc['sp_obl',names[i]] = sp_obl        
			sys_df.loc['sp_prc',names[i]] = sp_prc
			sys_df.loc['sp_lon',names[i]] = 0
		elif int(runprops.get('dynamicstoincludeflags')[i]) == 2:
			j2 = paramdf['j2r2_'+str(i+1)][0]
			c22 = paramdf['c22r2_'+str(i+1)][0]
			sp_rate = paramdf['sprate_'+str(i+1)][0]
			sp_obl = paramdf['spinc_'+str(i+1)][0]
			sp_prc = paramdf['splan_'+str(i+1)][0]
			sp_lon = paramdf['spaop_'+str(i+1)][0]      
			sys_df.loc['axis',names[i]] = axis      
			sys_df.loc['j2r2',names[i]] = j2     
			sys_df.loc['c22r2',names[i]] = c22    
			sys_df.loc['sp_rate',names[i]] = sp_rate
			sys_df.loc['sp_obl',names[i]] = sp_obl
			sys_df.loc['sp_prc',names[i]] = sp_prc
			sys_df.loc['sp_lon',names[i]] = sp_lon          
		if i > 0:
			sma = paramdf['sma_'+str(i+1)][0]
			ecc = paramdf['ecc_'+str(i+1)][0]
			inc = paramdf['inc_'+str(i+1)][0]
			lan = paramdf['lan_'+str(i+1)][0]
			aop = paramdf['aop_'+str(i+1)][0]
			mea = paramdf['mea_'+str(i+1)][0]
			sys_df.loc['sma',names[i]] = sma
			sys_df.loc['ecc',names[i]] = ecc
			sys_df.loc['inc',names[i]] = inc
			sys_df.loc['lan',names[i]] = lan
			sys_df.loc['aop',names[i]] = aop
			sys_df.loc['mea',names[i]] = mea
    
	#t_arr = times
	N = len(sys_df.columns)
    
	j2_sum = sum(sys_df.loc["j2r2",:].values.flatten())
	names = list(sys_df.columns)
	t_arr = np.array([0])
	totaltime = 50*365*24*3600
	t_arr = np.arange(0,totaltime,3600*24)
   
    
	tol = runprops.get("spinny_tolerance")

	#print(sys_df)
	if N == 2 and j2_sum == 0.00 and runprops.get('includesun') == 0:
		print(paramdf)        
		kepler_system = kepler_2body(paramdf,t_arr,runprops)
		kepler_df = kepler_system[0]
		names = kepler_system[1]
		spinny_plot(kepler_df, names,runprops)
	elif runprops.get('includesun') == 0:
		system = build_spinny_ns(sys_df,runprops)
		spinny = evolve_spinny_ns(system[0],system[1],system[2],system[3],system[4],system[5],t_arr,tol,runprops)
		s_df = spinny[0]
		names = spinny[2]
		spinny_plot(s_df,names, runprops)
	else: 
		system = build_spinny(sys_df, runprops)
		spinny = evolve_spinny(system[0],system[1],system[2],system[3],system[4],system[5],t_arr,runprops)
		        
		s_df = spinny[0]
		#print(s_df)
		#s_df = pd.DataFrame()
		#s_df['X_Pos_'+]=
		#s_df['']=
		#s_df['']=
		#s_df['']=
		#s_df['']=
		#s_df['']=        
		names = spinny[1]
		spinny_plot(s_df,names, runprops)
    
	#Model_DeltaLong, Model_DeltaLat, fakeobsdf = mm_likelihood.mm_chisquare(paramdf, fakeobsdf, runprops, geo_obj_pos, gensynth = True)
	DeltaLong_Model, DeltaLat_Model, fakeobsdf = mm_likelihood.mm_chisquare(paramdf, fakeobsdf, runprops, geo_obj_pos, gensynth = True)

	modelx = np.empty((nobjects-1, fakeobsdf.shape[0]))
	modely = np.empty((nobjects-1, fakeobsdf.shape[0]))

	x = np.empty((nobjects-1, obsdf.shape[0]))
	xe = np.empty((nobjects-1, obsdf.shape[0]))
	y = np.empty((nobjects-1, obsdf.shape[0]))
	ye = np.empty((nobjects-1, obsdf.shape[0]))

	colorcycle = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']
	#markercycle = ["x","+"]

	name_dict = runprops.get("names_dict")
	objectnames = []
	for i in name_dict.values():
		objectnames.append(i)

	fig = plt.figure()
	for i in range(1,nobjects):
		modelx[i-1,:] = DeltaLat_Model[i-1]
		modely[i-1,:] = DeltaLong_Model[i-1]

		x[i-1,:] = obsdf["DeltaLat_" + objectnames[i]].values
		xe[i-1,:] = obsdf["DeltaLat_" + objectnames[i] + "_err"].values
		y[i-1,:] = obsdf["DeltaLong_" + objectnames[i]].values
		ye[i-1,:] = obsdf["DeltaLong_" + objectnames[i] + "_err"].values

		plt.plot(modelx[i-1,:], modely[i-1,:], color = colorcycle[i], label = objectnames[i], linewidth = 0.5, alpha = 0.5)
		plt.errorbar(x[i-1,:], y[i-1,:], xerr = xe[i-1,:], yerr = ye[i-1,:], fmt = "ko", ms = 2)

	plt.axis('equal')
	plt.xlabel("Delta Latitude")
	plt.ylabel("Delta Longitude")
	plt.legend()

	plt.savefig("best_astrometry.pdf", format = "pdf")
	plt.close()

	obspdf = PdfPages("observations.pdf")

	modelpos = [modelx,modely]
	objpos = [x,y]
	objposerr = [xe,ye]
	labels = ["dLat", "dLong"]

	for i in range(1,nobjects):
		for j in range(2):
			plt.figure()
			plt.errorbar(time_arr, objpos[j][i-1,:], yerr = objposerr[j][i-1,:], fmt = "ko", ms = 2)
			plt.plot(time_arr, modelpos[j][i-1,:], colorcycle[i], linewidth = 0.75, alpha = 0.75, label = objectnames[i])
			plt.xlabel("Time (SJD)")
			plt.ylabel(labels[j])
			plt.legend()
			obspdf.savefig()

	obspdf.close()
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

                
#Actually build the plots here
#====================================================================================================
import glob, os

if 'results' in os.getcwd():
    getData = ReadJson('runprops.txt')
else:
    getData = ReadJson('most_recent_runprops.txt')
runprops = getData.outProps()
objname = runprops.get("objectname")

if not 'results' in os.getcwd():
	os.chdir('../../../results/'+objname+'/')
	results = max(glob.glob(os.path.join(os.getcwd(), '*/')), key=os.path.getmtime)
	os.chdir(results)

backend = emcee.backends.HDFBackend('chain.h5')
    
fit_scale = pd.read_csv('fit_scale.csv',index_col=0)
float_names = runprops.get('float_names')
obsdf = pd.read_csv(objname+'_obs_df.csv',index_col=0)
geo_obj_pos = pd.read_csv('geocentric_'+objname+'_position.csv',index_col=0)
fixed_df = pd.read_csv('fixed_df.csv',index_col=0)
total_df_names = runprops.get('total_df_names')

plots(backend, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names)


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
