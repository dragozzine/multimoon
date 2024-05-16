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

import matplotlib.animation as animation


def update(data):
    points.set_ydata(data)
    return points

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data
import glob,os

#def create(objname,filename):
#if 'results' in filename:
        #os.chdir(filename)
getData = ReadJson('runprops.txt')
#else:
        #getData = ReadJson('most_recent_runprops.txt')
runprops = getData.outProps()
objname = runprops.get("objectname")

    #if not 'results' in os.getcwd():
    #    os.chdir('../../../results/'+objname+'/')
    #    results = max(glob.glob(os.path.join(os.getcwd(), '*/')), key=os.path.getmtime)
    #    os.chdir(results)

backend = emcee.backends.HDFBackend('chain.h5')

init = pd.read_csv(objname+'_init_guess.csv',index_col=0)
print(init)
    
    
best_llhoods = pd.read_csv('best_likelihoods.csv')
best = best_llhoods.iloc[-(1)]
print(best)
if runprops.get('numobjects') == 2:
    ind = ['mass_1','mass_2','sma_2','ecc_2','aop_2','inc_2','lan_2','mea_2']
elif runprops.get('numobjects') == 3:
    ind = ['mass_1','mass_2','sma_2','ecc_2','aop_2','inc_2','lan_2','mea_2','mass_3','sma_3','ecc_3','aop_3','inc_3','lan_3','mea_3']

cols = ['mean','stddev']
best_init = pd.DataFrame(columns = cols, index = ind)
best_init['mean']['mass_1'] = best['mass_1']
best_init['mean']['mass_2'] = best['mass_2']
best_init['mean']['sma_2'] = best['sma_2']
best_init['mean']['ecc_2'] = best['ecc_2']
best_init['mean']['aop_2'] = best['aop_2']
best_init['mean']['inc_2'] = best['inc_2']
best_init['mean']['lan_2'] = best['lan_2']
best_init['mean']['mea_2'] = best['mea_2']
best_init['mean']['mass_3'] = best['mass_3']
best_init['mean']['sma_3'] = best['sma_3']
best_init['mean']['ecc_3'] = best['ecc_3']
best_init['mean']['aop_3'] = best['aop_3']
best_init['mean']['inc_3'] = best['inc_3']
best_init['mean']['lan_3'] = best['lan_3']
best_init['mean']['mea_3'] = best['mea_3']
best_init['stddev']['mass_1'] = init['stddev']['mass_1']
best_init['stddev']['mass_2'] = init['stddev']['mass_2']
best_init['stddev']['sma_2'] = init['stddev']['sma_2']
best_init['stddev']['ecc_2'] = init['stddev']['ecc_2']
best_init['stddev']['aop_2'] = init['stddev']['aop_2']
best_init['stddev']['inc_2'] = init['stddev']['inc_2']
best_init['stddev']['lan_2'] = init['stddev']['lan_2']
best_init['stddev']['mea_2'] = init['stddev']['mea_2']
best_init['stddev']['mass_3'] = init['stddev']['mass_3']
best_init['stddev']['sma_3'] = init['stddev']['sma_3']
best_init['stddev']['ecc_3'] = init['stddev']['ecc_3']
best_init['stddev']['aop_3'] = init['stddev']['aop_3']
best_init['stddev']['inc_3'] = init['stddev']['inc_3']
best_init['stddev']['lan_3'] = init['stddev']['lan_3']
best_init['stddev']['mea_3'] = init['stddev']['mea_3']
    #print(runprops)
runs_dir = runprops.get('objectname')
best_init.to_csv(runs_dir+'_init_guess_from_best.csv')
#print(best_init)
