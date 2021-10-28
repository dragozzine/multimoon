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
from tqdm import tqdm
import matplotlib.colors as colors
import functools

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

#chain = (nwalkers, nlink, ndim)

def sample_deltas(i, draws, names, fixed_df, total_df_names, fit_scale, names_dict, runprops, nobj, obsdf, geo_obj_pos):
		paramdf = mm_param.from_fit_array_to_param_df(draws[i,:].flatten(), names, fixed_df, total_df_names, fit_scale, names_dict, runprops)[0]
		dlong = np.zeros((runprops.get('numobjects')-1, 19))
		dlat = np.zeros((runprops.get('numobjects')-1, 19))
		resid_long = np.zeros((runprops.get('numobjects')-1, 19))
		resid_lat = np.zeros((runprops.get('numobjects')-1, 19))
		#print(paramdf.iloc[:,:-nobj].values)		
		drawparams = paramdf.iloc[:,:-nobj].values
		#print(paramdf)
		DeltaLong_Model, DeltaLat_Model, obsdf = mm_likelihood.mm_chisquare(paramdf, obsdf, runprops, geo_obj_pos, gensynth = True)
		chisq_total, residuals = mm_likelihood.mm_chisquare(paramdf, obsdf, runprops, geo_obj_pos)
		#print('Residuals ',residuals)        
		for j in range(0,runprops.get('numobjects')-1):
			dlong[j,:] = DeltaLong_Model[j]
			dlat[j,:] = DeltaLat_Model[j]
			resid_long[j,:] = residuals[2*j]
			resid_lat[j,:] = residuals[2*j+1]
		#print("dlong",dlong)
		#print("dlat",dlat)
		return dlong, dlat, drawparams, resid_long, resid_lat       

def posterior(sampler, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names, pool):
	numdraws = 1000

	# Getting log likelihood posterior values and flatchain for use throughout
	burnin = int(runprops.get('nburnin'))
	clusterburn = int(runprops.get('clustering_burnin'))
	thin_plots = int(runprops.get('nthinning'))
	flatchain = sampler.get_chain(discard=int(burnin/thin_plots+clusterburn/thin_plots),flat = True, thin=thin_plots)
	print(flatchain.shape, 'shape')
	all_llhoods = sampler.get_log_prob(discard=int(burnin/thin_plots+clusterburn/thin_plots),flat = True, thin=thin_plots)
	#ind = np.argmax(llhoods)
	#params = flatchain[ind,:].flatten()

	# Getting parameter names
	names = []
	for i in float_names:
		names.append(i)
	names_dict = runprops.get("names_dict")

	# Choose random draws from the flatchain
	drawsindex = np.random.randint(flatchain.shape[0], size = numdraws)
	draws = flatchain[drawsindex,:]
	#all_llhoods = sampler.get_log_prob(discard=int(burnin+clusterburn),flat = True, thin=thin_plots)
	llhoods = all_llhoods[drawsindex]

	# Get time arrays
	converttimes = ["2021-10-01","2022-09-30"]
	t = Time(converttimes)
	timesdic = {'start': t.isot[0], 'stop': t.isot[1], 'step': '1d'}

	# Make a geocentric position file
	#geo_obj_pos = mm_make_geo_pos.mm_make_geo_pos(objname, timesdic, runprops, True)

	# Creating a fake observtions data frame
	times = geo_obj_pos.values[:,0].flatten()
	fakeobsdf = obsdf.loc[[0,1],:]
	for i in range(len(times)):
		if i == 0 or i == 1:
			fakeobsdf.iloc[i,0] = times[i]
		fakeobsdf = fakeobsdf.append(fakeobsdf.iloc[-1,:])
		fakeobsdf['time'].iloc[-1] = times[i]
	fakeobsdf = fakeobsdf.iloc[2:]

	# Creating arrays to hold outputs
	dlong = np.zeros((draws.shape[0], runprops.get('numobjects')-1, times.size))
	dlat = np.zeros((draws.shape[0], runprops.get('numobjects')-1, times.size))
	resid_long = np.zeros((draws.shape[0], runprops.get('numobjects')-1, times.size))
	resid_lat = np.zeros((draws.shape[0], runprops.get('numobjects')-1, times.size))

	# Holding paramvalues
	nobj = runprops.get('numobjects')
	print(mm_param.from_fit_array_to_param_df(draws[0,:].flatten(), names, fixed_df, total_df_names, fit_scale, names_dict, runprops)[0])   
	ndims  = mm_param.from_fit_array_to_param_df(draws[0,:].flatten(), names, fixed_df, total_df_names, fit_scale, names_dict, runprops)[0].iloc[:,:-nobj].size
	print(ndims)
	paramnames  = mm_param.from_fit_array_to_param_df(draws[0,:].flatten(), names, fixed_df, total_df_names, fit_scale, names_dict, runprops)[0].columns.tolist()[0:-nobj]
	print(paramnames)
	drawparams = np.zeros((ndims, numdraws))

	# Looping to get model values
	deltas = functools.partial(sample_deltas, draws=draws, names=names, fixed_df=fixed_df, total_df_names=total_df_names, fit_scale=fit_scale, names_dict=names_dict, runprops=runprops, nobj=nobj, obsdf=obsdf, geo_obj_pos=geo_obj_pos)
	x = tqdm(range(draws.shape[0]))
	data = pool.map(deltas, x)

	dlong = np.zeros((draws.shape[0],2,19))
	dlat = np.zeros((draws.shape[0],2,19))
	resid_long = np.zeros((draws.shape[0],2,19))
	resid_lat = np.zeros((draws.shape[0],2,19))
	#print("DATA:    ", data[i][3])
	#print(data[i][4])
	for i in range(len(data)):          
		dlong[i] = data[i][0]
		dlat[i] = data[i][1]
		drawparams[:,i] = data[i][2]
		resid_long[i] = data[i][3]
		resid_lat[i] = data[i][4]
        
	'''
	for i in tqdm(range(draws.shape[0])):
		paramdf = mm_param.from_fit_array_to_param_df(draws[i,:].flatten(), names, fixed_df, total_df_names, fit_scale, names_dict, runprops)[0]
		#print(paramdf.iloc[:,:-nobj].values)		
		drawparams[:,i] = paramdf.iloc[:,:-nobj].values
		#print(paramdf)
		DeltaLong_Model, DeltaLat_Model, obsdf = mm_likelihood.mm_chisquare(paramdf, obsdf, runprops, geo_obj_pos, gensynth = True)
		for j in range(1,runprops.get('numobjects')):
			dlong[i,j-1,:] = DeltaLong_Model[j-1]
			dlat[i,j-1,:] = DeltaLat_Model[j-1]
		#print("dlong",dlong)
		#print("dlat",dlat)
	'''
	# Now collapse the arrays with a std call
	dlongstd = np.std(dlong,axis = 0)
	dlatstd = np.std(dlat,axis = 0)
	dlongmean = np.mean(dlong,axis = 0)
	dlatmean = np.mean(dlat,axis = 0)
	print(dlongstd.shape)
	print(dlatstd.shape)


	totaldf = pd.DataFrame(drawparams.T, columns = paramnames)
	#print(totaldf)


	# Calculate average (mean for now) error in the real data
	name_dict = runprops.get("names_dict")
	objectnames = []
	for i in name_dict.values():
		objectnames.append(i)
	typicalerror = np.zeros((2,runprops.get('numobjects')-1))
	for i in range(1,runprops.get('numobjects')):
		typicalerror[0,i-1] = np.median(obsdf["DeltaLong_" + objectnames[i] + "_err"].values.flatten())
		typicalerror[1,i-1] = np.median(obsdf["DeltaLat_" + objectnames[i] + "_err"].values.flatten())

	# Plot dlong vs dlat with color for j2
	from matplotlib.backends.backend_pdf import PdfPages

	#print("deltas ", dlong[:,0,0], dlat[:,0,0],llhoods)
	#print(dlong.shape, len(dlong[0,0,:]))
#	exit()


	#===================================Here we create the residuals heat map========================================
	#for i in range(len(dlong[0,0,:])):
	#	chisquare_total, residuals = mm_likelihood.mm_chisquare(drawparams[i,:], obsdf, runprops, geo_obj_pos)
	#	resid_list[i] = residuals
	residpdf = PdfPages("resid_map.pdf")
	for i in range(1, runprops.get('numobjects')-1):        
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
		#print(nobjects, np.array(residuals).shape, objectnames)  

		print('plotting ', i, ' ',objectnames[i])
		plt.hist2d(resid_long[:,1,i], resid_lat[:,1,i], bins=40, range=[[-4.0, 4.0], [-3.0, 3.0]],label = objectnames[i], edgecolors = None)
		plt.xlabel("Delta Longitude")
		plt.ylabel("Delta Latitude")
		plt.axis("equal")
		plt.legend()
	residpdf.close()        

	#==============================================We create the posterior.pdf====================================================
	predictionspdf = PdfPages("posterior.pdf")
	markers = ["o","D","^"]    
	for i in range(len(dlong[0,0,:])):
		plt.figure()
		plt.axis("equal")        
		plt.scatter(0,0, color = "black")
		for j in range(runprops.get('numobjects')-1):        
			plt.scatter(dlong[:,j,i], dlat[:,j,i], c=llhoods, cmap = "coolwarm",marker=markers[j])
			plt.errorbar(np.median(dlong[:,j,i]), np.median(dlat[:,j,i]), xerr = typicalerror[0,j], yerr = typicalerror[1,j], ecolor = "red")
		plt.xlabel("dLon")
		plt.ylabel("dLat")
		#plt.xlim(-0.5, 0.5)
		#plt.ylim(-0.5, 0.5)
		plt.title("JD "+str(obsdf['time'][i]))
		color_bar = plt.colorbar()
		color_bar.set_alpha(1)
		color_bar.draw_all()
		color_bar.set_label('Log-Likelihood')
		predictionspdf.savefig()
	predictionspdf.close()


	#==============================================We create the brightness.pdf====================================================
	if runprops.get("photo_offset"):
		mass_rat = totaldf['mass_2']/totaldf['mass_1']
		bright_rat = abs(totaldf['f_val_1'])*mass_rat**(2/3)
		hubble_sep_arc = 2.1*10**5*5.5*10**(-7)/2.4  
		brightnesspdf = PdfPages("brightness.pdf")
		for i in range(len(dlong[0,0,:])):
			plt.figure()
			plt.axis('equal')
			plt.scatter(0,0, color = "black")
			plt.scatter(dlong[:,0,i], dlat[:,0,i], c=bright_rat, cmap = "coolwarm", norm=colors.LogNorm())
			plt.xlabel("dLon")
			plt.ylabel("dLat")
			#plt.xlim(-0.2, 0.2)
			#plt.ylim(-0.2, 0.2)
			plt.title("JD "+str(obsdf['time'][i]))
			color_bar = plt.colorbar()
			color_bar.set_alpha(1)
			color_bar.draw_all()
			color_bar.set_label('Brightness ratio')
			brightnesspdf.savefig()
		brightnesspdf.close()

		#==========================================We create the brightness_seperation.pdf================================================= 
		brightnessratpdf = PdfPages("brightness_seperation.pdf")
		for i in range(len(dlong[0,0,:])):
			plt.figure()   
			plt.axis('equal')
			sep = np.sqrt(dlong[:,0,i]**2+dlat[:,0,i]**2)        
			plt.scatter(sep, bright_rat, c=llhoods, cmap = "coolwarm")
			plt.axvline(x=hubble_sep_arc, color='r')
			plt.axvline(x=0.2, color='b')
			plt.xlabel("total separation")
			plt.ylabel("brightness ratio")   
			#plt.xlim(0, 0.3)
			#plt.ylim(0, 0.15)
			plt.title("JD "+str(obsdf['time'][i]))
			color_bar = plt.colorbar()
			color_bar.set_label('Log-Likelihood')
			color_bar.set_alpha(1)
			color_bar.draw_all()
			brightnessratpdf.savefig()
		brightnessratpdf.close()
            
            
#Actually build the plots here
#====================================================================================================
import glob, os

if __name__ == '__main__':

    from schwimmbad import MPIPool
    with MPIPool() as pool:
    
        if not pool.is_master():
            pool.wait()
            sys.exit(0)

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

        posterior(backend, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names, pool)
