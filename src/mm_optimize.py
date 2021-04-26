#
#	mm_optimize.py
#
#	Uses an optimizer to get a better initial guess
#
#	Benjamin Proudfoot
#	09/28/20
#

import numpy as np
import pandas as pd
import random
import commentjson as json
import emcee
import scipy
from tqdm import tqdm
import functools
from datetime import datetime

import mm_likelihood

def neg_log_prob(float_params, float_names, fixed_df, total_df_names, fit_scale, runprops, obsdf, geo_obj_pos, best_llhoods):
	out =  -1*mm_likelihood.log_probability(float_params, float_names, fixed_df, total_df_names, fit_scale, runprops, obsdf, geo_obj_pos, best_llhoods)
	#print('neg_log_prob: ',out)
	return out

def op_map(i, p0, float_names, fixed_df, total_df_names, fit_scale, runprops, obsdf, geo_obj_pos, best_llhoods):
    return scipy.optimize.minimize(neg_log_prob, p0[i,:], args = (float_names,fixed_df.iloc[[i]],total_df_names, fit_scale, runprops, obsdf, geo_obj_pos, best_llhoods), method = runprops.get('opt_method'), tol=runprops.get("spinny_tolerance"),options={'maxiter': runprops.get('opt_iter')})

def mm_optimize(nwalkers, p0, float_names, fixed_df, total_df_names, fit_scale, runprops, obsdf, geo_obj_pos, best_llhoods, pool):
	#print(p0)    
	optimization = functools.partial(op_map, p0=p0, float_names=float_names, fixed_df = fixed_df, total_df_names=total_df_names, fit_scale=fit_scale, runprops=runprops, obsdf=obsdf, geo_obj_pos=geo_obj_pos, best_llhoods=best_llhoods)
	x = tqdm(range(nwalkers))
	begin = datetime.now()    
	data = pool.map(optimization, x)
	print(datetime.now()-begin)    
	for i in range(len(data)):
		p0[i,:]= data[i].x
    
	#for i in tqdm(range(nwalkers)):
	#	res = scipy.optimize.minimize(neg_log_prob, p0[i,:], args = (float_names,fixed_df.iloc[[i]],total_df_names, fit_scale, runprops, obsdf, geo_obj_pos, best_llhoods), method = "Powell")
	#	p0[i,:] = res.x
	#	if res.success:
	#		print("Optimization successful")
	#		print("result: ", res.x)
	#		print()
	#	else:
	#		print("Optimization unsuccessful")
	return p0
