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
    
getData = ReadJson('runprops.txt')

runprops = getData.outProps()
objname = runprops.get("objectname")

backend = emcee.backends.HDFBackend('chain.h5')

chain = backend.get_chain(flat = False, thin=thin_plots)  

pd.to_csv('chain_thinned.csv')