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

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data
    
    
if 'results' in os.getcwd():
    getData = ReadJson('runprops.txt')
#else:
    #getData = ReadJson('most_recent_runprops.txt')
    runprops = getData.outProps()
    objname = runprops.get("objectname")

#if not 'results' in os.getcwd():
#    objname = runprops.get('objectname')
#    os.chdir('../../../results/'+objname+'/')
#    results = max(glob.glob(os.path.join(os.getcwd(), '*/')), key=os.path.getmtime)
#    os.chdir(results)
    

    backend = emcee.backends.HDFBackend('chain.h5')

    fit_scale = pd.read_csv('fit_scale.csv',index_col=0)
    float_names = runprops.get('float_names')
    obsdf = pd.read_csv(objname+'_obs_df.csv',index_col=0)
    geo_obj_pos = pd.read_csv('geocentric_'+objname+'_position.csv',index_col=0)
    fixed_df = pd.read_csv('fixed_df.csv',index_col=0)
    total_df_names = runprops.get('total_df_names')
    
    thin_plots = runprops.get('thin_plots')
    
    chain = sampler.get_chain(flat = False, thin=thin_plots) 
    
    np.savetxt('thin_chain.txt',chain)