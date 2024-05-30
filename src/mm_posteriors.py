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


class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

#chain = (nwalkers, nlink, ndim)
   

def posterior(sampler, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names):
	numdraws = 100

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

	# Holding paramvalues
	nobj = runprops.get('numobjects')
	print(mm_param.from_fit_array_to_param_df(draws[0,:].flatten(), names, fixed_df, total_df_names, fit_scale, names_dict, runprops)[0])   
	ndims  = mm_param.from_fit_array_to_param_df(draws[0,:].flatten(), names, fixed_df, total_df_names, fit_scale, names_dict, runprops)[0].iloc[:,:-nobj].size
	print(ndims)
	paramnames  = mm_param.from_fit_array_to_param_df(draws[0,:].flatten(), names, fixed_df, total_df_names, fit_scale, names_dict, runprops)[0].columns.tolist()[0:-nobj]
	print(paramnames)
	drawparams = np.zeros((ndims, numdraws))

	# Looping to get model values
    
	
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
	predictionspdf = PdfPages("posterior.pdf")
	for i in range(len(dlong[0,0,:])):
		plt.figure()
		plt.scatter(0,0, color = "black")
		plt.scatter(dlong[:,0,i], dlat[:,0,i], c=llhoods, cmap = "coolwarm")
		plt.errorbar(np.median(dlong[:,0,i]), np.median(dlat[:,0,i]), xerr = typicalerror[0,0], yerr = typicalerror[1,0], ecolor = "red")
        
		plt.scatter(dlong[:,1,i], dlat[:,1,i], c=llhoods, cmap = "coolwarm",marker="D")
		plt.errorbar(np.median(dlong[:,1,i]), np.median(dlat[:,1,i]), xerr = typicalerror[0,1], yerr = typicalerror[1,1], ecolor = "red")
		plt.xlabel("dLon")
		plt.ylabel("dLat")
		plt.xlim(-0.5, 0.5)
		plt.ylim(-0.5, 0.5)
		plt.title("JD "+str(obsdf['time'][i]))
		color_bar = plt.colorbar()
		color_bar.set_alpha(1)
		color_bar.draw_all()
		predictionspdf.savefig()
	predictionspdf.close()

	if runprops.get("photo_offset"):
		brightnesspdf = PdfPages("brightness.pdf")
		for i in range(len(dlong[0,0,:])):
			plt.figure()
			plt.scatter(0,0, color = "black")
			plt.scatter(dlong[:,0,i], dlat[:,0,i], c=totaldf['f_val_1'], cmap = "coolwarm", norm=colors.LogNorm())
			plt.xlabel("dLon")
			plt.ylabel("dLat")
			plt.xlim(-0.2, 0.2)
			plt.ylim(-0.2, 0.2)
			plt.title("JD "+str(obsdf['time'][i]))
			color_bar = plt.colorbar()
			color_bar.set_alpha(1)
			color_bar.draw_all()
			brightnesspdf.savefig()
		brightnesspdf.close()
        
		mass_rat = totaldf['mass_2']/totaldf['mass_1']
		bright_rat = totaldf['f_val_1']*mass_rat**(2/3)
		hubble_sep_arc = 2.1*10**5*5.5*10**(-7)/2.4  
		brightnessratpdf = PdfPages("brightness_ratio.pdf")
		for i in range(len(dlong[0,0,:])):
			plt.figure()            
			sep = np.sqrt(dlong[:,0,i]**2+dlat[:,0,i]**2)        
			plt.scatter(sep, bright_rat, c=llhoods, cmap = "coolwarm")
			plt.axvline(x=hubble_sep_arc, color='r')
			plt.axvline(x=0.2, color='b')
			plt.xlabel("total separation")
			plt.ylabel("brightness ratio")   
			plt.xlim(0, 0.3)
			plt.ylim(0, 0.15)
			plt.title("JD "+str(obsdf['time'][i]))
			color_bar = plt.colorbar()
			color_bar.set_alpha(1)
			color_bar.draw_all()
			brightnessratpdf.savefig()
		brightnessratpdf.close()
            
            
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

posterior(backend, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names)