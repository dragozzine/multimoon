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

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

#chain = (nwalkers, nlink, ndim)

def posteriordraw(sampler, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names):
    numdraws = 100      # How many posterior draws to take
    years = 5        # How many years to run the integrations for
    tstep = 1          # Time step size in days

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

    # Get time arrays and set constants
    epoch = runprops.get("epoch_SJD")
    converttimes = [epoch,epoch+(365.25*years)]
    t = Time(converttimes, format = "jd")
    timesdic = {'start': t.isot[0], 'stop': t.isot[1], 'step': str(tstep)+'d'}

    # Make a geocentric position file
    geo_obj_pos = mm_make_geo_pos.mm_make_geo_pos(objname, timesdic, runprops, True)
    
    # Creating a fake observtions data frame
    times = geo_obj_pos.values[:,0].flatten()
    times_earth = geo_obj_pos.values[:,1].flatten()
    fakeobsdf = obsdf.loc[[0,1],:]
    for i in range(len(times)):
        if i == 0 or i == 1:
            fakeobsdf.iloc[i,0] = times[i]
        fakeobsdf = fakeobsdf.append(fakeobsdf.iloc[-1,:])
        fakeobsdf['time'].iloc[-1] = times[i]
    fakeobsdf = fakeobsdf.iloc[2:]

    # Holding paramvalues
    nobj = runprops.get('numobjects')
    ndims  = mm_param.from_fit_array_to_param_df(draws[0,:].flatten(), names, fixed_df, total_df_names, fit_scale, names_dict, runprops)[0].iloc[:,:-nobj].size
    paramnames  = mm_param.from_fit_array_to_param_df(draws[0,:].flatten(), names, fixed_df, total_df_names, fit_scale, names_dict, runprops)[0].columns.tolist()[0:-nobj]
    drawparams = np.zeros((ndims, numdraws))
    
    # Making an empty dataframe to enable concatenation down the line
    outdf = pd.DataFrame()

    # Looping to get model values
    for i in tqdm(range(draws.shape[0])):
        paramdf = mm_param.from_fit_array_to_param_df(draws[i,:].flatten(), names, fixed_df, total_df_names, fit_scale, names_dict, runprops)[0]
        drawparams[:,i] = paramdf.iloc[:,:-nobj].values
        s_df = mm_likelihood.mm_chisquare(paramdf, fakeobsdf, runprops, geo_obj_pos, dfout = True)
        s_df["Times"] = times_earth
        #print(s_df)
        outdf = pd.concat([outdf, s_df], axis = 0)
        #print(outdf)
        
    
    # Now save all the arrays
    outdf.to_csv("draw_integrations.csv")
    pd.DataFrame(drawparams.T, columns = paramnames).to_csv("draw_params.csv")

#Actually run the code here
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

posteriordraw(backend, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names)

# Now I'll give a short example for how to extract the data, it extracts the mean anomaly

# First read in the two dataframes that hold the integrations and the parameters
integrations = pd.read_csv("draw_integrations.csv", index_col = 0)
params = pd.read_csv("draw_params.csv", index_col = 0)

# Find the number of draws and number of steps taken
ndraws = params.shape[0]
nsteps = integrations.idxmax()["Times"]

# Start making a plot
plt.figure()

# get the time array. only need to do it once since they're all the same
t = integrations["Times"].iloc[:(nsteps+1)].to_numpy()
tp = Time(t, format = "jd")

# Now get the mean anomaly out. It goes [draw1/time1, draw2/time1, ... drawN/timeM, drawN+1/timeM ...]
# You need to do some slick array slicing after the fact to get it out nicely
mea = integrations["ecc_Weywot"].loc[np.arange(nsteps+1)]

# Loop over the individual posterior draws to make a plot
for i in range(ndraws):
    # Now time for the array slicing
    mea_draw = mea[i::ndraws]
    plt.plot(tp.datetime,mea_draw, color = "black", alpha = 0.1)

plt.xlabel("Time")
plt.ylabel("Mean anomaly")
plt.savefig("testout.png")











