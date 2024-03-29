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
import functools

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

#chain = (nwalkers, nlink, ndim)

def sample_deltas(i, draws, names, fixed_df, total_df_names, fit_scale, names_dict, runprops, nobj, fakeobsdf, geo_obj_pos, dlong, dlat):
		paramdf = mm_param.from_fit_array_to_param_df(draws[i,:].flatten(), names, fixed_df, total_df_names, fit_scale, names_dict, runprops)[0]
       
		drawparams = paramdf.iloc[:,:-nobj].values
		DeltaLong_Model, DeltaLat_Model, fakeobsdf = mm_likelihood.mm_chisquare(paramdf, fakeobsdf, runprops, geo_obj_pos, gensynth = True)
		length = len(DeltaLong_Model)        
		dlong = np.zeros((runprops.get('numobjects')-1, length))
		dlat = np.zeros((runprops.get('numobjects')-1, length))
		for j in range(1,runprops.get('numobjects')):
			dlong[j-1,:] = DeltaLong_Model[j-1]
			dlat[j-1,:] = DeltaLat_Model[j-1]
		return dlong, dlat, drawparams        

def predictions(sampler, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names, pool):
	numdraws = 1000

	# Getting log likelihood posterior values and flatchain for use throughout
	burnin = int(runprops.get('nburnin'))
	clusterburn = int(runprops.get('clustering_burnin'))
	thin_plots = int(runprops.get('nthinning'))
	flatchain = sampler.get_chain(discard=int(burnin/thin_plots+clusterburn/thin_plots),flat = True, thin=thin_plots)
	print(flatchain.shape, 'shape')
	llhoods = sampler.get_log_prob(discard=int(burnin/thin_plots+clusterburn/thin_plots),flat = True, thin=thin_plots)

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
	llhoods = llhoods[drawsindex]

	# Get time arrays
	converttimes = ["2023-10-01","2024-09-30"]
	t = Time(converttimes)
	timesdic = {'start': t.isot[0], 'stop': t.isot[1], 'step': '6h'}

	# Make a geocentric position file
	geo_obj_pos = mm_make_geo_pos.mm_make_geo_pos(objname, timesdic, runprops, True)

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

	deltas = functools.partial(sample_deltas, draws=draws, names=names, fixed_df=fixed_df, total_df_names=total_df_names, fit_scale=fit_scale, names_dict=names_dict, runprops=runprops, nobj=nobj, fakeobsdf=fakeobsdf, geo_obj_pos=geo_obj_pos, dlong=dlong, dlat=dlat)
	x = tqdm(range(draws.shape[0]))
	data = pool.map(deltas, x)
	length = len(data)    

	dlong = np.zeros((draws.shape[0],2,length))
	dlat = np.zeros((draws.shape[0],2,length))
	for i in range(len(data)):          
		dlong[i] = data[i][0]
		dlat[i] = data[i][1]
		drawparams[:,i] = data[i][2]
        
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

	# Now create info gain arrays
	infogain = np.zeros((runprops.get('numobjects')-1, times.size))
	infogain2 = np.zeros((runprops.get('numobjects')-1, times.size))
	g_gain = np.zeros((runprops.get('numobjects')-1,2))   
	g_gains_ind = np.zeros((runprops.get('numobjects')-1,10))
	for i in range(1,runprops.get('numobjects')):
		infogain[i-1,:] = np.sqrt( (dlongstd[i-1,:]/typicalerror[0,i-1])**2 + (dlatstd[i-1,:]/typicalerror[1,i-1])**2 )
		g_gain[i-1,0] = times[np.argmax(infogain[i-1,:])]
		g_gain[i-1,1] = np.argmax(infogain[i-1,:])
		g_gains_ind[i-1,:] = (-infogain[i-1,:]).argsort()[:10]        


	colorcycle = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']
	fig = plt.figure(figsize = (12.8,4.8))
	t = Time(times, format = "jd")
	mat_times = matplotlib.dates.datestr2num(t.iso, default=None)
	for i in range(1,runprops.get('numobjects')):
		plt.plot_date(mat_times, infogain[i-1,:].flatten(), "-", color = colorcycle[i-1], label = objectnames[i], alpha = 0.5)
		plt.title("Dates of greatest gain: JD "+str(g_gain[:,0]))        

	plt.xlabel("Time")
	plt.ylabel("Info gained")
	plt.legend()
	plt.savefig("predictions.pdf", format = "pdf")
	plt.close()

	print('dlatmean',dlatmean)
	inds = g_gains_ind[0].astype(int) 
	print('dlatmean',dlatmean[0,inds])
	gains = pd.DataFrame(columns=['times','dates','Hercules_Delta_Lat','Hercules_Delta_Long','Hercules_angle_tot','Dysnomia_Delta_Lat','Dysnomia_Delta_Long','Dysnomia_angle_tot','infogain_val'])
	gains['times'] = times[inds]
	gains['dates'] = t.iso[inds]

	# DS TODO: remove eris specific code

	gains['Hercules_Delta_Lat'] = dlatmean[0,inds]
	gains['Hercules_Delta_Long'] = dlongmean[0,inds]
	gains['Hercules_angle_tot'] = np.sqrt(dlatmean[0,inds]**2+dlongmean[0,inds]**2)
	gains['Dysnomia_Delta_Lat'] = dlatmean[1,inds]
	gains['Dysnomia_Delta_Long'] = dlongmean[1,inds]
	gains['Dysnomia_angle_tot'] = np.sqrt(dlatmean[1,inds]**2+dlongmean[1,inds]**2)
	gains['infogain_val'] = infogain[1,inds]
	print(gains)
	gains.to_csv('predictions_df.csv')    
	# Plot dlong vs dlat with color for j2
	from matplotlib.backends.backend_pdf import PdfPages

	predictionspdf = PdfPages("predictions_params.pdf")
	print("Making predictions params")
	for i in range(len(paramnames)):
		plt.figure()
		plt.axis("equal")
		plt.scatter(0,0, color = "black")
		plt.scatter(dlong[:,0,int(g_gain[1,1])], dlat[:,0,int(g_gain[1,1])], c = totaldf[paramnames[i]], edgecolor = None, alpha = 0.5, s = 10, cmap = "coolwarm")
		plt.errorbar(np.median(dlong[:,0,int(g_gain[1,1])]), np.median(dlat[:,0,int(g_gain[1,1])]), xerr = typicalerror[0,0], yerr = typicalerror[1,0], ecolor = "red")
		plt.scatter(dlong[:,1,int(g_gain[1,1])], dlat[:,1,int(g_gain[1,1])], c = totaldf[paramnames[i]], edgecolor = None, alpha = 0.5, s = 10, cmap = "coolwarm",marker='D')
		plt.errorbar(np.median(dlong[:,1,int(g_gain[1,1])]), np.median(dlat[:,1,int(g_gain[1,1])]), xerr = typicalerror[0,1], yerr = typicalerror[1,1], ecolor = "red")
		plt.xlabel("dLon")
		plt.ylabel("dLat")
		plt.title(paramnames[i])
		color_bar = plt.colorbar()
		color_bar.set_alpha(1)
		color_bar.draw_all()
		color_bar.set_label(paramnames[i])
		predictionspdf.savefig()
	predictionspdf.close()

    
	pospdf = PdfPages("posterior_prediction.pdf")
	plt.figure()
	plt.axis("equal")        
	plt.scatter(0,0, color = "black")
       
	plt.scatter(dlong[:,0,int(g_gain[1,1])], dlat[:,0,int(g_gain[1,1])], c=llhoods, cmap = "coolwarm")
	plt.errorbar(np.median(dlong[:,0,int(g_gain[1,1])]), np.median(dlat[:,0,int(g_gain[1,1])]), xerr = typicalerror[0,0], yerr = typicalerror[1,0], ecolor = "red")
       
	plt.scatter(dlong[:,1,int(g_gain[1,1])], dlat[:,1,int(g_gain[1,1])], c=llhoods, cmap = "coolwarm",marker="D")
	plt.errorbar(np.median(dlong[:,1,int(g_gain[1,1])]), np.median(dlat[:,1,int(g_gain[1,1])]), xerr = typicalerror[0,1], yerr = typicalerror[1,1], ecolor = "red")
	plt.xlabel("dLon")
	plt.ylabel("dLat")
	plt.title("JD: "+str(g_gain[1,0]))
	color_bar = plt.colorbar()
	color_bar.set_alpha(1)
	color_bar.draw_all()
	color_bar.set_label('Log-Likelihood')
	pospdf.savefig()
	pospdf.close()

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

        predictions(backend, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names, pool)
