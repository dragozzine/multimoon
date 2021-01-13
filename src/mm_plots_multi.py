import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import corner
import numpy as np
import pandas as pd
import emcee
import sys
import mm_likelihood
from astropy.time import Time
import commentjson as json
import mm_param
import mm_make_geo_pos

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

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
	undo_spin = np.zeros(runprops.get('numobjects')-1)
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
    
	full_chain = sampler.get_chain(flat = False, thin=100)
	flatchain = sampler.get_chain(discard=(burnin+clusterburn),flat = True, thin=100)
	print(flatchain.shape, full_chain.shape)    
	fit = []

	for i in fit_scale.columns:
		name = i
		if type(name) != str:
			name = name[0]
		if name in float_names:
			val = fit_scale.loc[0, i]
			fit.append(val)
                  
	chain = sampler.get_chain(discard=(burnin+clusterburn),flat = False, thin=100)
	print(chain.shape)
	numparams = chain.shape[2]
	numwalkers = chain.shape[1]
	numgens = chain.shape[0]
 
	#First fit the flatchain with the fit total_df_names  
	#print(numwalkers, numgens, numparams)
	#print(flatchain, len(flatchain),len(flatchain[0]), fit)
	fchain = np.zeros((numgens*numwalkers,numparams))    
	for i in range(numgens*numwalkers):
		row = np.zeros(numparams)        
		for j in range(numparams):
			#print(i, ',', j)
			val = flatchain[i][j]*fit[j]
			row[j] = val
		#print(row)
		if runprops.get('transform'):
			for b in range(runprops.get('numobjects')-1):           
				if undo_ecc_aop[b]:
					aop_new = row[int(ecc_aop_index[b*2+1])]
					ecc_new = row[int(ecc_aop_index[b*2])]
					pomega = np.arctan2(ecc_new,aop_new)*180/np.pi
					if pomega < 0:
						pomega = pomega%360
					row[int(ecc_aop_index[b*2+1])] = pomega
					row[int(ecc_aop_index[b*2])] = ecc_new/np.sin(pomega*np.pi/180)
                               
				if undo_inc_lan[b]:
					inc_new = row[int(inc_lan_index[b*2])]
					lan_new = row[int(inc_lan_index[b*2+1])]
					lan = np.arctan2(inc_new,lan_new)*180/np.pi
					if lan < 0:
						lan = lan%360
					row[int(inc_lan_index[b*2+1])] = lan
					inc = np.arctan2(inc_new,np.sin(lan*np.pi/180))*2*180/np.pi
					if inc < 0:
						inc = inc%180
					row[int(inc_lan_index[b*2])] = inc
                    
				if undo_spin[b]:
					spinc_new = row[int(spin_index[b*2])]
					splan_new = row[int(spin_index[b*2+1])]
					splan = np.arctan2(spinc_new,splan_new)*180/np.pi
					if splan < 0:
						splan = splan%360
					row[int(spin_index[b*2+1])] = splan
					spinc = np.arctan2(spinc_new,np.sin(splan*np.pi/180))*2*180/np.pi
					if spinc < 0:
						spinc = spinc%180
					row[int(spin_index[b*2])] = spinc
                
				if undo_lambda[b]:
					mea_new = row[int(lambda_index[b*2])]
					pomega = row[int(lambda_index[b*2+1])]
                
					mea = mea_new-pomega
					if mea < 0:
						mea = mea%360
					elif mea > 360:
						mea = mea%360
					row[int(lambda_index[b*2])] = mea
                    
				if undo_pomega[b]:
					lan = row[int(pomega_index[b*2+1])]
					pomega = row[int(pomega_index[b*2])]
					aop = pomega-lan
					if aop < 0:
						aop = aop%360
					elif mea > 360:
						aop = aop%360
					row[int(pomega_index[b*2])] = aop
                    
			if undo_masses[0]:
				mass_1 = row[int(masses_index[0])]
				mass_2 = row[int(masses_index[1])]
				row[int(masses_index[1])] = mass_2-mass_1
			elif undo_masses[1]:
				mass_1 = row[int(masses_index[0])]
				mass_2 = row[int(masses_index[1])]
				mass_3 = row[int(masses_index[2])]
				row[int(masses_index[2])] = mass_3-mass_2 
				row[int(masses_index[1])] = mass_2-mass_1 
		fchain[i] = np.array(row)

	flatchain = np.array(fchain)
	#print(flatchain)
    
	names = []
	for i in float_names:
		names.append(i)
	transform_names = np.copy(names)        

	#Now fit the chain
	extra_chain = []
	cchain = np.zeros((numgens,numwalkers, numparams))    
	for i in range(numgens):
		for j in range(numwalkers):
			row = []
			for k in range(numparams):
				val = chain[i][j][k]*fit[k]
				cchain[i][j][k] = val
			if runprops.get('transform'):
				for b in range(runprops.get('numobjects')-1):
					if undo_ecc_aop[b]:
						np.append(transform_names,['ecc*sin(pomega)'+str(b+1),'ecc*cos(pomega)'+str(b+1)])
						aop_new = cchain[i][j][int(ecc_aop_index[b*2+1])]
						ecc_new = cchain[i][j][int(ecc_aop_index[b*2])]
						extra_chain.append(ecc_new)
						extra_chain.append(aop_new)
						pomega = np.arctan2(ecc_new,aop_new)*180/np.pi
						if pomega < 0:
							pomega = pomega%360
						cchain[i][j][int(ecc_aop_index[b*2+1])] = pomega
						cchain[i][j][int(ecc_aop_index[b*2])] = ecc_new/np.sin(pomega/180*np.pi)
					if undo_inc_lan[b]:    
						np.append(transform_names,['tan(inc/2)*sin(lan)_'+str(b+1),'tan(inc/2)*cos(lan)_'+str(b+1)])   
						inc_new = cchain[i][j][int(inc_lan_index[b*2])]
						lan_new = cchain[i][j][int(inc_lan_index[b*2+1])]
						extra_chain.append(inc_new)
						extra_chain.append(lan_new)
						lan = np.arctan2(inc_new,lan_new)*180/np.pi
						if lan < 0:
							lan = lan%360
						cchain[i][j][int(inc_lan_index[b*2+1])] = lan
						inc = np.arctan2(inc_new,np.sin(lan*np.pi/180))*2*180/np.pi
						if inc < 0:
							inc = inc%180
						cchain[i][j][int(inc_lan_index[b*2])] = inc
					if undo_spin[b]:
						np.append(transform_names,['tan(spinc/2)*sin(splan)_'+str(b+1),'tan(spinc/2)*cos(splan)_'+str(b+1)])   
						spinc_new = cchain[i][j][int(spin_index[b*2])]
						splan_new = cchain[i][j][int(spin_index[b*2+1])]
						extra_chain.append(spinc_new)
						extra_chain.append(splan_new)
						splan = np.arctan2(spinc_new,splan_new)*180/np.pi
						if splan < 0:
							splan = splan%360
						cchain[i][j][int(spin_index[b*2+1])] = lan
						spinc = np.arctan2(spinc_new,np.sin(splan*np.pi/180))*2*180/np.pi
						if spinc < 0:
							spinc = spinc%180
						cchain[i][j][int(spin_index[b*2])] = spinc
					if undo_lambda[b]:
						np.append(transform_names,['lambda(pomega+mea)_'+str(b+1)])   
						mea_new = cchain[i][j][int(lambda_index[b*2])]
						pomega = cchain[i][j][int(lambda_index[b*2+1])]
						extra_chain.append(mea_new)
						mea = mea_new-pomega
						if mea < 0:
							mea = mea%360
						if mea > 360:
							mea = mea%360
						cchain[i][j][int(lambda_index[b*2])] = mea
					if undo_pomega[b]:
						np.append(transform_names,['pomega_'+str(b+1)])   
						lan = cchain[i][j][int(pomega_index[b*2+1])]
						pomega = cchain[i][j][int(pomega_index[b*2])]
						extra_chain.append(pomega)
						aop = pomega-lan
						if aop < 0:
							aop = aop%360
						if aop > 360:
							aop = aop%360                        
						cchain[i][j][int(pomega_index[b*2])] = pomega-lan  
				if undo_masses[0]:
					np.append(transform_names,['mass1+mass2'])   
					mass_1 = cchain[i][j][int(masses_index[0])]
					mass_2 = cchain[i][j][int(masses_index[1])]
					extra_chain.append(mass_2)
					cchain[i][j][int(masses_index[1])] = mass_2-mass_1
				elif undo_masses[1]:
					np.append(transform_names,['mass1+mass2+mass3'])   
					mass_1 = cchain[i][j][int(masses_index[0])]
					mass_2 = cchain[i][j][int(masses_index[1])]
					mass_3 = cchain[i][j][int(masses_index[2])]
					extra_chain.append(mass_3)
					cchain[i][j][int(masses_index[2])] = mass_3-mass_2 
					cchain[i][j][int(masses_index[1])] = mass_2-mass_1 
                        
	cchain = np.array(cchain)

	oldchain = chain
	chain = cchain

	# Make corner plot
	#plt.rc('text', usetex=True)
	fig = 0
	#if runprops.get("objectname") == "testcases" and runprops.get("unseenmoon"):
	#	getData = ReadJson("../runs/"+objname+"/"+runprops.get('run_file')+"/runprops_gensynth.txt")
	#	synthrunprops = getData.outProps()
	#	truths = []
	#
	#	params_dict = synthrunprops.get("params_dict")
	#
	#	for k in params_dict.values():
	#		truths.append(k)
	#
	#	fig = corner.corner(flatchain, labels = names, bins = 40, show_titles = True, 
	#			    plot_datapoints = False, color = "blue", fill_contours = True,
	#			    title_fmt = ".4e", truths = truths)
	print(flatchain)
	fig = corner.corner(flatchain, labels = names, bins = 40, show_titles = True, 
			    plot_datapoints = False, color = "blue", fill_contours = True,
			    title_fmt = ".4e")
	fig.tight_layout(pad = 1.08, h_pad = 0, w_pad = 0)
	for ax in fig.get_axes():
		ax.tick_params(axis = "both", labelsize = 20, pad = 0.5)
	fname = "corner.pdf"       
	fig.savefig(fname, format = 'pdf')
	plt.close("all")
	#plt.rc('text', usetex=False)
	

	# Now make the walker plots
	from matplotlib.backends.backend_pdf import PdfPages

	walkerpdf = PdfPages("walkers.pdf")

	for i in range(numparams):
		plt.figure()
		for j in range(numwalkers):
			plt.plot(np.reshape(chain[0:numgens,j,i], numgens))
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
	fullgens = int(numgens/100+burnin/100+clusterburn/100)
	#print(fullgens)
	for i in range(numparams):
		plt.figure()
		for j in range(numwalkers):
			plt.plot(np.reshape(full_chain[0:fullgens,j,i], fullgens))
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
	old_fchain = sampler.get_chain(flat=True)
	llhoods = sampler.get_log_prob(discard=(burnin+clusterburn),flat = True, thin=100)
	sigsdf = pd.DataFrame(columns = ['-3sigma','-2sigma','-1sigma','median','1sigma','2sigma','3sigma', 'mean'], index = transform_names)
	j = 0
	for i in range(len(flatchain[0])):
		num = flatchain[:,i]
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
	filename = 'sigsdf.csv'    
	sigsdf.to_csv(filename)
    
	# Likelihood plots    
	likelihoodspdf = PdfPages("likelihoods.pdf")

	for i in range(numparams):
		plt.figure(figsize = (9,9))
		plt.subplot(221)
		plt.hist(flatchain[:,i].flatten(), bins = 40, histtype = "step", color = "black")
		plt.subplot(223)
		plt.scatter(flatchain[:,i].flatten(), llhoods.flatten(),
			    c = np.mod(np.linspace(0,llhoods.size - 1, llhoods.size), numwalkers),
			    cmap = "nipy_spectral", edgecolors = "none", rasterized=True)
		plt.xlabel(names[i])
		plt.ylabel("Log(L)")
		plt.subplot(224)
		plt.hist(llhoods.flatten(), bins = 40, orientation = "horizontal", 
			 histtype = "step", color = "black")
		likelihoodspdf.attach_note(names[i])
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
	paramdf = mm_param.from_fit_array_to_param_df(paraminput, paramnames, fixed_df, total_df_names, fit_scale, names_dict, runprops)

	#paramdf = pd.DataFrame(paraminput).transpose()
	#paramdf.columns = paramnames


	#print(paramdf)
#Currently this function call sends an error in the case of leaving any necessary value floating, since paramdf will be incomplete 
	chisquare_total, residuals = mm_likelihood.mm_chisquare(paramdf, obsdf, runprops, geo_obj_pos)
	#print(chisquare_total, residuals)

	colorcycle = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']

	xvals = np.linspace(-1.0,1.0,num=1000)
	circle = np.sqrt(1 - xvals**2)

	plt.figure()
	plt.plot(xvals, circle, color = "black")
	plt.plot(xvals,-circle, color = "black")
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
			plt.plot(times, modelpos[j][i-1,:], colorcycle[i], linewidth = 0.75, alpha = 0.75, label = objectnames[i])
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
