#
#	unpack.py
#
#	Unpacks chain.h5 files
#
#	Benjamin Proudfoot
#	12/2/21
#

import numpy as np
import pandas as pd
import os
from tqdm import tqdm
from scipy import interpolate
import sys
import commentjson as json
import emcee
class ReadJson(object):
    def __init__(self, filename):
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

def unpack(resultdir, thin_plot = 1):

	# Loading in all necessary files for unpacking chain
	getData = ReadJson(resultdir+'runprops.txt')
	runprops = getData.outProps()
	objname = runprops.get("objectname")
	chain = emcee.backends.HDFBackend(resultdir+'chain.h5')
	fit_scale = pd.read_csv(resultdir+'fit_scale.csv',index_col=0)
	float_names = runprops.get('float_names')
	obsdf = pd.read_csv(resultdir+objname+'_obs_df.csv',index_col=0)
	geo_obj_pos = pd.read_csv(resultdir+'geocentric_'+objname+'_position.csv',index_col=0)
	fixed_df = pd.read_csv(resultdir+'fixed_df.csv',index_col=0)
	total_df_names = runprops.get('total_df_names')

	# Begin to unpack chain
	nthinning = runprops.get("nthinning")
	undo_ecc_aop = np.zeros(runprops.get('numobjects')-1)
	undo_ecc_aop[:] = False
	ecc_aop_index = np.zeros((runprops.get('numobjects')-1)*2)
	undo_inc_lan = np.zeros(runprops.get('numobjects')-1)
	undo_inc_lan[:] = False
	undo_spin = np.zeros(runprops.get('numobjects'))
	undo_spin[:] = False
	spin_index = np.zeros((runprops.get('numobjects'))*2)
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
			undo_spin[i+1] = True
			spin_index[2*(i+1)] = float_names.index('spinc_'+str(i+2))
			spin_index[2*(i+1)+1] = float_names.index('splan_'+str(i+2))
		if 'mea_'+str(i+2) in float_names and 'aop_'+str(i+2) in float_names:
			undo_lambda[i] = True
			lambda_index[2*i] = float_names.index('mea_'+str(i+2))
			lambda_index[2*i+1] = float_names.index('aop_'+str(i+2))
		if 'aop_'+str(i+2) in float_names and 'lan_'+str(i+2) in float_names:
			undo_pomega[i] = True
			pomega_index[2*i] = float_names.index('aop_'+str(i+2))
			pomega_index[2*i+1] = float_names.index('lan_'+str(i+2))
	if 'spinc_1' in float_names and 'splan_1' in float_names:
		undo_spin[0] = True
		spin_index[0] = float_names.index('spinc_1')
		spin_index[1] = float_names.index('splan_1')
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

	if runprops.get('thin_run'):
		burnin = int(runprops.get('nburnin')/runprops.get('nthinning'))
		clusterburn = int(runprops.get('clustering_burnin')/runprops.get('nthinning'))
	else:
		burnin = int(runprops.get('nburnin'))
		clusterburn = int(runprops.get('clustering_burnin'))

	#thin_plots = runprops.get('thin_plots') 
#	chain = sampler.get_chain(discard=int(burnin+clusterburn),flat = False, thin=thin_plots)  
	chain = sampler.get_chain(flat = False, thin=thin_plots)  
	fit = []

	for i in fit_scale.columns:
		name = i
		if type(name) != str:
			name = name[0]
		if name in float_names:
			val = fit_scale.loc[0, i]
			fit.append(val)

	# Getting final values for the shape of the chain
	#shortchain = sampler.get_chain(discard=int(burnin+clusterburn),flat = False, thin=thin_plots)
	shortchain = sampler.get_chain(flat = False, thin=thin_plots)
	numparams = shortchain.shape[2]
	numwalkers = shortchain.shape[1]
	numgens = shortchain.shape[0]
	del shortchain
 
	# Take chain "fit" values and make them into real values
	for i in range(numparams):
		chain[:,:,i] = chain[:,:,i]*fit[i]

    
	fitparam_chain = np.zeros((1,numwalkers,numgens))

	print(fitparam_chain.shape)    
	fitparam_names = []    
	# Now de-transform the chain
	print("Starting un transformations")
	if runprops.get("transform"):
		for b in range(runprops.get('numobjects')-1):
			if undo_ecc_aop[b]:
				aop_new = chain[:,:,int(ecc_aop_index[b*2+1])]
				ecc_new = chain[:,:,int(ecc_aop_index[b*2])]
				#print(aop_new.T.shape, np.array([aop_new.T]).shape)
				fitparam_chain = np.concatenate((fitparam_chain, np.array([aop_new.T])),axis=0)
				fitparam_chain = np.concatenate((fitparam_chain, np.array([ecc_new.T])),axis=0)                
				fitparam_names.append('equinoctial_k_'+str(b+1))
				fitparam_names.append('equinoctial_h_'+str(b+1))
				pomega = (np.arctan2(ecc_new,aop_new)*180/np.pi)%360
				chain[:,:,int(ecc_aop_index[b*2+1])] = pomega
				chain[:,:,int(ecc_aop_index[b*2])] = ecc_new/np.sin(pomega/180*np.pi)
			if undo_inc_lan[b]:
				inc_new = chain[:,:,int(inc_lan_index[b*2])]
				lan_new = chain[:,:,int(inc_lan_index[b*2+1])]
				fitparam_chain = np.concatenate((fitparam_chain, np.array([inc_new.T])),axis=0)
				fitparam_chain = np.concatenate((fitparam_chain, np.array([lan_new.T])),axis=0)
				fitparam_names.append('equinoctial_q_'+str(b+1))
				fitparam_names.append('equinoctial_p_'+str(b+1))
				lan = (np.arctan2(inc_new,lan_new)*180/np.pi)%360
				chain[:,:,int(inc_lan_index[b*2+1])] = lan
				inc = (np.arctan2(inc_new,np.sin(lan*np.pi/180))*2*180/np.pi)%180
				chain[:,:,int(inc_lan_index[b*2])] = inc
			if undo_lambda[b]:
				mea_new = chain[:,:,int(lambda_index[b*2])]
				pomega = chain[:,:,int(lambda_index[b*2+1])]
				fitparam_chain = np.concatenate((fitparam_chain, np.array([mea_new.T])),axis=0)
				fitparam_chain = np.concatenate((fitparam_chain, np.array([pomega.T])),axis=0)
				fitparam_names.append('lambda_'+str(b+1))
				fitparam_names.append('pomega_'+str(b+1))
				mea = (mea_new-pomega)%360
				chain[:,:,int(lambda_index[b*2])] = mea
			if undo_pomega[b]:
				lan = chain[:,:,int(pomega_index[b*2+1])]
				pomega = chain[:,:,int(pomega_index[b*2])]
				aop = (pomega-lan)%360
				chain[:,:,int(pomega_index[b*2])] = aop
		for b in range(runprops.get('numobjects')):
			if undo_spin[b]:
				spinc_new = chain[:,:,int(spin_index[b*2])]
				splan_new = chain[:,:,int(spin_index[b*2+1])]
				fitparam_chain = np.concatenate((fitparam_chain, np.array([spinc_new.T])),axis=0)
				fitparam_chain = np.concatenate((fitparam_chain, np.array([splan_new.T])),axis=0)
				fitparam_names.append('spin_equinoctial_p_'+str(b+1))
				fitparam_names.append('spin_equinoctial_q_'+str(b+1))
				splan = (np.arctan2(spinc_new,splan_new)*180/np.pi)%360
				chain[:,:,int(spin_index[b*2+1])] = splan
				#spinc = (np.arctan2(spinc_new,np.sin(splan*np.pi/180))*2*180/np.pi)%180
				#print('spinc_new ',spinc_new)
				#print('splan_new ',splan_new)
				spinc = (np.arccos(spinc_new/np.sin(splan*np.pi/180))*2*180/np.pi)%180
				#print('spinc ',spinc)
				#print('splan ',splan)
				chain[:,:,int(spin_index[b*2])] = spinc
		if undo_masses[0]:
			mass_1 = chain[:,:,int(masses_index[0])]
			mass_2 = chain[:,:,int(masses_index[1])]
			fitparam_chain = np.concatenate((fitparam_chain, np.array([mass_2.T])),axis=0)
			fitparam_names.append('mass1+2')
			chain[:,:,int(masses_index[1])] = (mass_2-mass_1)/(10**18)
			chain[:,:,int(masses_index[0])] = (mass_1)/(10**18)  
			#print('hallo')
		elif undo_masses[1]:
			mass_1 = chain[:,:,int(masses_index[0])]
			mass_2 = chain[:,:,int(masses_index[1])]
			mass_3 = chain[:,:,int(masses_index[2])]
			#print(mass_1,mass_2, mass_3)            
			fitparam_chain = np.concatenate((fitparam_chain, np.array([mass_2.T])),axis=0)
			fitparam_chain = np.concatenate((fitparam_chain, np.array([mass_3.T])),axis=0)
			fitparam_names.append('mass1+2')
			fitparam_names.append('mass1+2+3')
			chain[:,:,int(masses_index[2])] = (mass_3-mass_2)/10**18
			chain[:,:,int(masses_index[1])] = (mass_2-mass_1)/10**18
			chain[:,:,int(masses_index[0])] = (mass_1)/10**18

	fitparam_chain = np.delete(fitparam_chain,0,0)
	fitparam_chain = fitparam_chain.T   
	fitparam_chain = fitparam_chain[int(burnin+clusterburn + thin_plots - 1) :: 1]
        
	print("Un transforming done")

	# Cutting up chain
	full_chain = np.copy(chain)
	chain = chain[int(burnin+clusterburn + thin_plots - 1) :: 1]
	print(chain.shape)

	# Flattening the chain based on method in emcee
	s = list(chain.shape[1:])
	s[0] = np.prod(chain.shape[:2])
	s2 = list(fitparam_chain.shape[1:])
	s2[0] = np.prod(fitparam_chain.shape[:2])
	#print(masses_index,s2, fitparam_chain.shape)    
	flatchain = chain.reshape(s)
	fitparam_chain = fitparam_chain.reshape(s2)    
	#print('flatchain and fitparam chain shape\n',flatchain.shape, fitparam_chain.shape)

	#print('flatchain[:,0] and flatchain[:,1]\n',flatchain[:,0],'\n',flatchain[:,1],'\n',flatchain[:,2],'\n',flatchain[:,3],'\n',flatchain[:,4],'\n',flatchain[:,5],'\n',flatchain[:,6])
	# Getting parameter names
	names = []
	#print(float_names)    
	for i in float_names:
		names.append(i)

	return flatchain, names

