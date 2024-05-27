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
import os.path

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

#chain = (nwalkers, nlink, ndim)

def posteriordraw(sampler, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names):
    numdraws = 2       # How many posterior draws to take
    years = 1          # How many years to run the integrations for
    tstep = 1          # Time step size in days must be an integer
    
    # Check if a past posterior draw has been completed
    if not os.path.isfile("draw_params.csv"):
        # Now draw from the posterior if the file doesn't exist
        print("Drawing from the posterior")
        
        # Getting flatchain
        burnin = int(runprops.get('nburnin'))
        clusterburn = int(runprops.get('clustering_burnin'))
        thin_plots = int(runprops.get('nthinning'))
        flatchain = sampler.get_chain(discard=int(burnin/thin_plots+clusterburn/thin_plots),flat = True, thin=thin_plots)

        # Choose random draws from the flatchain
        drawsindex = np.random.randint(flatchain.shape[0], size = numdraws)
        draws = flatchain[drawsindex,:]
    # If a draw has already been taken and saved, this just loads it in and uses it
    else:
        print("Using existing posterior draw")
        drawdf = pd.read_csv("draw_params.csv", index_col = 0)
        numdraws = drawdf.shape[0]
        
    # Getting parameter names, this is required for MultiMoon to run
    names = []
    for i in float_names:
        names.append(i)
    names_dict = runprops.get("names_dict")

    # Get time arrays set up, this will be inputs for the Horizons query in future steps
    # Integrations will always run from the chosen epoch (in runprops)
    epoch = runprops.get("epoch_SJD")
    converttimes = [epoch,epoch+(365.25*years)]
    t = Time(converttimes, format = "jd")
    timesdic = {'start': t.isot[0], 'stop': t.isot[1], 'step': str(tstep)+'d'}

    # Make a geocentric position file, using MultiMoon's infrastructure
    geo_obj_pos = mm_make_geo_pos.mm_make_geo_pos(objname, timesdic, runprops, True)
    
    # Creating a fake observtions data frame
    # This acts as a dummy input to satisfy input requirements
    times = geo_obj_pos.values[:,0].flatten()
    times_earth = geo_obj_pos.values[:,1].flatten()
    fakeobsdf = obsdf.loc[[0,1],:]
    for i in range(len(times)):
        if i == 0 or i == 1:
            fakeobsdf.iloc[i,0] = times[i]
        fakeobsdf = fakeobsdf.append(fakeobsdf.iloc[-1,:])
        fakeobsdf['time'].iloc[-1] = times[i]
    fakeobsdf = fakeobsdf.iloc[2:]

    # Prepping to store each paramdf
    # draw params holds the actual paramdf of each posterior draw
    nobj = runprops.get('numobjects')
    drawparams = pd.DataFrame([], columns = ["mass_1", "mass_2","sma_2","ecc_2","aop_2","inc_2","lan_2","mea_2",
                                             "j2r2_1","c22r2_1","spinc_1","splan_1","spaop_1","sprate_1",
                                             "j2r2_2","c22r2_2","spinc_2","splan_2","spaop_2","sprate_2",
                                             "pbad","logjitter","name_1","name_2"])
    
    # Making an empty dataframe to enable concatenation down the line
    # This will store all the integrations
    outdf = pd.DataFrame()
    
    # Make a dummy paramdf with all possible parameters
    # This acts as a way to add any parameters to the paramdf that aren't there originally
    dummy = pd.DataFrame([], columns = ["mass_1", "mass_2","sma_2","ecc_2","aop_2","inc_2","lan_2","mea_2",
                                        "j2r2_1","c22r2_1","spinc_1","splan_1","spaop_1","sprate_1",
                                        "j2r2_2","c22r2_2","spinc_2","splan_2","spaop_2","sprate_2",
                                        "pbad","logjitter","name_1","name_2"])

    # Looping to get model values
    for i in tqdm(range(numdraws)):
        # Get paramdf from running the appropriate command within multimoon
        # Only used when drawing from the posterior for the first time
        if not os.path.isfile("draw_params.csv"):
            paramdf = mm_param.from_fit_array_to_param_df(draws[i,:].flatten(), names, fixed_df, total_df_names, 
                                                          fit_scale, names_dict, runprops)[0]
        # If using a previous draw, it just uses the stored version.
        # This allows the draws to be altered between runs in excel!
        # For example, you can add J2 for the secondary to explore its spin dynamics
        # Warning! Changing values may make the fit to the data bad!!!
        else:
            paramdf = drawdf.iloc[[i]]
        
        # Do some dataframe bookkeeping here
        paramdf = pd.concat([dummy, paramdf]).fillna(value = 0.0)
        drawparams = pd.concat([drawparams, paramdf])
        
        # Actually do the integration here
        s_df = mm_likelihood.mm_chisquare(paramdf, fakeobsdf, runprops, geo_obj_pos, dfout = True)
        
        # Replace the times in s_df with observation time in days
        # s_df outputs time in seconds on a clock co-moving with the TNB
        s_df["Times"] = times_earth

        # Now concatenate the current integration with any past integrations
        # They stack on top of each other!
        outdf = pd.concat([outdf, s_df], axis = 0)
        
    # Now save all the arrays
    outdf.to_csv("draw_integrations.csv")
    if not os.path.isfile("draw_params.csv"):
        drawparams.to_csv("draw_params.csv")
    sys.exit()

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
# These have to be treated appropriately down the line to extract values correctly
integrations = pd.read_csv("draw_integrations.csv", index_col = 0)
params = pd.read_csv("draw_params.csv", index_col = 0)

# Find the number of draws and number of time steps taken
ndraws = params.shape[0]
nsteps = integrations.idxmax()["Times"]

# Start making a plot
plt.figure()

# get the time array. only need to do it once since they're all the same
t = integrations["Times"].iloc[:(nsteps+1)].to_numpy()
tp = Time(t, format = "jd")

# Now get the mean anomaly out. Integrations goes [draw1/time1, draw2/time1, ... drawN/timeM, drawN+1/timeM ...]
# You need to do some slick array slicing after the fact to get it out nicely
mea = integrations["ecc_Weywot"].loc[np.arange(nsteps+1)]

# Loop over the individual posterior draws to make a plot
for i in range(ndraws):
    # Now time for the array slicing.
    # The i::ndraws means start at i, then draw every ndraws indices to the end of the array
    mea_draw = mea[i::ndraws]
    plt.plot(tp.datetime, mea_draw, color = "black", alpha = 0.1)

plt.xlabel("Time")
plt.ylabel("Mean anomaly")
plt.savefig("testout.png")

