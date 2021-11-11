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
from mm_SPINNY.spinny_plots import astro_plot_time
from mm_SPINNY.spinny_generate import *
from mm_SPINNY.spinny_vector import *
from mm_SPINNY.spinny_nosun import *
from mm_SPINNY.mm_vpython import *
from mm_SPINNY.keplerian import *

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data
def save(s_df,names):

    save_yn = "N"
    if save_yn=="Y":
        print("Generating .csv...")
        t_current = ctime().replace(" ","_")
        filename = names[1]+"_SPINNY_"+t_current+".csv"
        s_df.to_csv("output/"+filename)
        print("SPINNY data saved to the output file as "+filename)
        plot_q(s_df, names)
    elif save_yn == "N":
        print("")
        plot_q(s_df, names)
    else:
        print("")
        print('Invalid Response.')
        return save(s_df,names)   
    
    
#chain = (nwalkers, nlink, ndim)

def plots(sampler, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names):
    # Here total_df_names is whatever file/object will have the run params
    objname = runprops.get("objectname")
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
    #print(os.getcwd())    
    
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

    thin_plots = runprops.get('thin_plots') 
#    chain = sampler.get_chain(discard=int(burnin+clusterburn),flat = False, thin=thin_plots)  
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
                spinc = (np.arctan2(spinc_new,np.sin(splan*np.pi/180))*2*180/np.pi)%180
                #spinc = (np.arccos(spinc_new/np.sin(splan*np.pi/180))*2*180/np.pi)%180
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


    # Getting parameter names
    names = []
    #print(float_names)    
    for i in float_names:
        names.append(i)


    # Getting log likelihood posterior values for use throughout
    all_llhoods = sampler.get_log_prob(discard=int(burnin+clusterburn),flat = True, thin=thin_plots)
    #print('llhoods shape',llhoods.shape)    
    ind = np.argmax(all_llhoods)
    params = flatchain[ind,:].flatten()
    from matplotlib.backends.backend_pdf import PdfPages
    #llhoods = sampler.get_log_prob(discard=int(burnin+clusterburn),flat = True, thin=thin_plots)
    flatchain = sampler.get_chain(flat = True)
    nobjects = runprops.get('numobjects')
    llhoods = sampler.get_log_prob(flat = True)
    ind = np.argmax(llhoods)
    params = flatchain[ind,:].flatten()
    
    
    
    #Get draws from the llhood and flatchain
    #===========================================================================================================
    numdraws = 4
    drawsindex = np.random.randint(flatchain.shape[0], size = numdraws)
    draws = flatchain[drawsindex,:]
    #all_llhoods = sampler.get_log_prob(discard=int(burnin+clusterburn),flat = True, thin=thin_plots)
    llhoods = llhoods[drawsindex]
    num = 0
    name_dict = runprops.get("names_dict")
    paramnames = names
    
    #Run this loop 20 times, which is the number of posterior draws we pull.
    sys_dfs = []
    for i in draws:
        paraminput = []
        num = num+1
        params = i.flatten()
        print('params ', params)
        for j in params:
            paraminput.append(j)    

        objectnames = []
        for j in range(runprops.get('numobjects')):
            objectnames.append(name_dict.get('name_'+str(j+1)))
        
        for values,keys in name_dict.items():
            for j in range(runprops.get('numobjects')):
                if str(j+1) in keys or 'offset' in keys: 
                    paraminput.append(values)
                    paramnames.append(keys)
    
        print('paraminput ',paraminput)
        #print(paramnames)
        names_dict = runprops.get("names_dict")    
        paramdf,fit_params = mm_param.from_fit_array_to_param_df(paraminput, paramnames, fixed_df, total_df_names, fit_scale, names_dict, runprops)
    
        print(paramdf)
        #print(draws)
    #Currently this function call sends an error in the case of leaving any necessary value floating, since paramdf will be incomplete 
        #chisquare_total, residuals = mm_likelihood.mm_chisquare(paramdf, obsdf, runprops, geo_obj_pos)
        # Astrometry plots
        time_arr = obsdf['time'].values.flatten()
        tmin = time_arr.min()
        tmax = time_arr.max()
    
        converttimes = [tmin,tmax]
        t = Time(converttimes, format = 'jd')
    
        timesdic = {'start': t.isot[0], 'stop': t.isot[1], 'step': '6h'}
        
        #geo_obj_pos = mm_make_geo_pos.mm_make_geo_pos(objname, timesdic, runprops, True)
        geo_obj_pos = pd.read_csv('geocentric_'+objname+'_position_analysis.csv')
    
        times = geo_obj_pos.values[:,0].flatten()
    
        fakeobsdf = obsdf.loc[[0,1],:]
        for j in range(len(times)):
            if j == 0 or j == 1:
                fakeobsdf.iloc[j,0] = times[j]
                # change row number?
            fakeobsdf = fakeobsdf.append(fakeobsdf.iloc[-1,:])
            fakeobsdf['time'].iloc[-1] = times[j]
        fakeobsdf = fakeobsdf.iloc[2:]
    
        names = list(runprops.get('names_dict').values())
        vals = ['mass','axis','j2r2','c22r2','sp_rate','sp_obl','sp_prc','sp_lon','sma','ecc','inc','lan','aop','mea']    
        sys_df = pd.DataFrame(columns=names,index=vals)
        sys_df[names[0]]['sma'] = 0
        sys_df[names[0]]['ecc'] = 0
        sys_df[names[0]]['inc'] = 0
        sys_df[names[0]]['lan'] = 0
        sys_df[names[0]]['aop'] = 0
        sys_df[names[0]]['mea'] = 0
    
        if runprops.get('includesun'):
            sys_df.loc['mass','Sun'] = paramdf['mass_0'][0]
            axis = list(runprops.get('axes_size').values())[0]
            j2 = 0
            c22 = 0
            sp_rate = 0
            sp_obl = 0
            sp_prc = 0
            sp_lon = 0      
            sys_df.loc['axis','Sun'] = axis      
            sys_df.loc['j2r2','Sun'] = j2     
            sys_df.loc['c22r2','Sun'] = c22    
            sys_df.loc['sp_rate','Sun'] = sp_rate
            sys_df.loc['sp_obl','Sun'] = sp_obl
            sys_df.loc['sp_prc','Sun'] = sp_prc
            sys_df.loc['sp_lon','Sun'] = sp_lon          
            sma = paramdf['sma_0'][0]
            ecc = paramdf['ecc_0'][0]
            inc = paramdf['inc_0'][0]
            lan = paramdf['lan_0'][0]
            aop = paramdf['aop_0'][0]
            mea = paramdf['mea_0'][0]
            sys_df.loc['sma','Sun'] = sma
            sys_df.loc['ecc','Sun'] = ecc
            sys_df.loc['inc','Sun'] = inc
            sys_df.loc['lan','Sun'] = lan
            sys_df.loc['aop','Sun'] = aop
            sys_df.loc['mea','Sun'] = mea
        
        for j in range(runprops.get('numobjects')):
            sys_df.loc['mass',names[j]] = paramdf['mass_'+str(j+1)][0]
            axis = list(runprops.get('axes_size').values())[j]
            if int(runprops.get('dynamicstoincludeflags')[j]) == 0:        
                sys_df.loc['axis',names[j]] = 0        
                sys_df.loc['j2r2',names[j]] = 0        
                sys_df.loc['c22r2',names[j]] = 0        
                sys_df.loc['sp_rate',names[j]] = 0        
                sys_df.loc['sp_obl',names[j]] = 0        
                sys_df.loc['sp_prc',names[j]] = 0        
                sys_df.loc['sp_lon',names[j]] = 0
            elif int(runprops.get('dynamicstoincludeflags')[j]) == 1:   
                j2 = paramdf['j2r2_'+str(j+1)][0]
                c22 = paramdf['c22r2_'+str(j+1)][0]
                sp_rate = paramdf['sprate_'+str(j+1)][0]
                sp_obl = paramdf['spinc_'+str(j+1)][0]
                sp_prc = paramdf['splan_'+str(j+1)][0]
                sp_lon = paramdf['spaop_'+str(j+1)][0]
                sys_df.loc['axis',names[j]] = axis
                sys_df.loc['j2r2',names[j]] = j2        
                sys_df.loc['c22r2',names[j]] = 0
                sys_df.loc['sp_rate',names[j]] = sp_rate
                sys_df.loc['sp_obl',names[j]] = sp_obl        
                sys_df.loc['sp_prc',names[j]] = sp_prc
                sys_df.loc['sp_lon',names[j]] = 0
            elif int(runprops.get('dynamicstoincludeflags')[j]) == 2:
                j2 = paramdf['j2r2_'+str(j+1)][0]
                c22 = paramdf['c22r2_'+str(j+1)][0]
                sp_rate = paramdf['sprate_'+str(j+1)][0]
                sp_obl = paramdf['spinc_'+str(j+1)][0]
                sp_prc = paramdf['splan_'+str(j+1)][0]
                sp_lon = paramdf['spaop_'+str(j+1)][0]      
                sys_df.loc['axis',names[j]] = axis      
                sys_df.loc['j2r2',names[j]] = j2     
                sys_df.loc['c22r2',names[j]] = c22    
                sys_df.loc['sp_rate',names[j]] = sp_rate
                sys_df.loc['sp_obl',names[j]] = sp_obl
                sys_df.loc['sp_prc',names[j]] = sp_prc
                sys_df.loc['sp_lon',names[j]] = sp_lon          
            if j > 0:
                sma = paramdf['sma_'+str(j+1)][0]
                ecc = paramdf['ecc_'+str(j+1)][0]
                inc = paramdf['inc_'+str(j+1)][0]
                lan = paramdf['lan_'+str(j+1)][0]
                aop = paramdf['aop_'+str(j+1)][0]
                mea = paramdf['mea_'+str(j+1)][0]
                sys_df.loc['sma',names[j]] = sma
                sys_df.loc['ecc',names[j]] = ecc
                sys_df.loc['inc',names[j]] = inc*np.pi/180
                sys_df.loc['lan',names[j]] = lan*np.pi/180
                sys_df.loc['aop',names[j]] = aop*np.pi/180
                sys_df.loc['mea',names[j]] = mea*np.pi/180
        
        #t_arr = times
        
        N = len(sys_df.columns)
        
        j2_sum = sum(sys_df.loc["j2r2",:].values.flatten())
        names = list(sys_df.columns)
        sys_dfs.append(sys_df)
        
    t_arr = np.array([0])
    totaltime = 50*365*24*3600
    t_arr = np.arange(0,totaltime,3600*24)
       
        
    tol = runprops.get("spinny_tolerance")
    
        #print(sys_df)
    runprops['animate'] = True
    runprops['animate_num'] = str(num)
    s_dfs = []
    nameses = []
    if N == 2 and j2_sum == 0.00 and runprops.get('includesun') == 0:
        print(paramdf)        
        kepler_system = kepler_2body(paramdf,t_arr,runprops)
        kepler_df = kepler_system[0]
        names = kepler_system[1]
        spinny_plot(kepler_df, names,runprops)
    elif runprops.get('includesun') == 0 and N==3:
        for i in sys_dfs:
            system = build_spinny_ns(i,runprops)
            spinny = evolve_spinny_ns(system[0],system[1],system[2],system[3],system[4],system[5],t_arr,tol,runprops)
            s_df = spinny[0]
            names = spinny[2]
            s_dfs.append(s_df)
            nameses.append(names)
        #astro_plot_time(s_df,names, runprops)
        spinny_plot_multiple(s_dfs,names,t_arr)
        
        #spinny_plot(s_df,names, runprops)
    elif runprops.get('includesun') == 1: 
        system = build_spinny(sys_dfs[0], runprops)
        spinny = evolve_spinny(system[0],system[1],system[2],system[3],system[4],system[5],t_arr,runprops)
        s_df = spinny[0]     
        names = spinny[1]
        spinny_plot(s_df,names, runprops)
     
    #Model_DeltaLong, Model_DeltaLat, fakeobsdf = mm_likelihood.mm_chisquare(paramdf, fakeobsdf, runprops, geo_obj_pos, gensynth = True)
        
    return()
    
#HEre we plot the time delayed astrometry with spinny
def plot_astro(sampler, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names):
    objname = runprops.get("objectname")
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
    #print(os.getcwd())    
    
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

    thin_plots = runprops.get('thin_plots') 
#    chain = sampler.get_chain(discard=int(burnin+clusterburn),flat = False, thin=thin_plots)  
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
                spinc = (np.arctan2(spinc_new,np.sin(splan*np.pi/180))*2*180/np.pi)%180
                #spinc = (np.arccos(spinc_new/np.sin(splan*np.pi/180))*2*180/np.pi)%180
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


    # Getting parameter names
    names = []
    #print(float_names)    
    for i in float_names:
        names.append(i)


    # Getting log likelihood posterior values for use throughout
    from matplotlib.backends.backend_pdf import PdfPages
    #llhoods = sampler.get_log_prob(discard=int(burnin+clusterburn),flat = True, thin=thin_plots)
    flatchain = sampler.get_chain(flat = True)
    nobjects = runprops.get('numobjects')
    llhoods = sampler.get_log_prob(flat = True)
    ind = np.argmax(llhoods)
    params = flatchain[ind,:].flatten()
    
    
    
    #Get draws from the llhood and flatchain

    num = 0
    name_dict = runprops.get("names_dict")
    paramnames = names
    
    paraminput = []
    num = num+1
    #params = i.flatten()
    #print('params ', params)
    for j in params:
        paraminput.append(j)    

    objectnames = []
    for j in range(runprops.get('numobjects')):
        objectnames.append(name_dict.get('name_'+str(j+1)))
        
    for values,keys in name_dict.items():
        for j in range(runprops.get('numobjects')):
            if str(j+1) in keys or 'offset' in keys: 
                paraminput.append(values)
                paramnames.append(keys)
    
    
    #print(paramnames)
    names_dict = runprops.get("names_dict")    
    paramdf,fit_params = mm_param.from_fit_array_to_param_df(paraminput, paramnames, fixed_df, total_df_names, fit_scale, names_dict, runprops)
    
    
    #print(draws)
    #Currently this function call sends an error in the case of leaving any necessary value floating, since paramdf will be incomplete 
    #chisquare_total, residuals = mm_likelihood.mm_chisquare(paramdf, obsdf, runprops, geo_obj_pos)
        # Astrometry plots
    time_arr = obsdf['time'].values.flatten()
    tmin = time_arr.min()
    tmax = time_arr.max()
    
    converttimes = [tmin,tmax]
    t = Time(converttimes, format = 'jd')
    
    timesdic = {'start': t.isot[0], 'stop': t.isot[1], 'step': '3h'}
        
    #geo_obj_pos = mm_make_geo_pos.mm_make_geo_pos(objname, timesdic, runprops, True)
    geo_obj_pos = pd.read_csv('geocentric_'+objname+'_position_analysis.csv')[-410:]
    #print(geo_obj_pos)
    #times = geo_obj_pos.values[:,0].flatten().copy()
    times = geo_obj_pos['kboTIME'].values.flatten().copy()
    
    fakeobsdf = obsdf.loc[[0,1],:]
    for j in range(len(times)):
        if j == 0 or j == 1:
            fakeobsdf.iloc[j,0] = times[j]
            # change row number?
        fakeobsdf = fakeobsdf.append(fakeobsdf.iloc[-1,:])
        #print('times[j] ', times[j])
        fakeobsdf['time'].iloc[-1] = times[j]
    fakeobsdf = fakeobsdf.iloc[2:]
    #astro_plot_time(param_df, names, runprops)
    
    
    DeltaLong_Model, DeltaLat_Model, fakeobsdf = mm_likelihood.mm_chisquare(paramdf, fakeobsdf, runprops, geo_obj_pos, gensynth = True)
    
    modelx = np.empty((nobjects-1, fakeobsdf.shape[0]))
    modely = np.empty((nobjects-1, fakeobsdf.shape[0]))
    
    x = np.empty((nobjects-1, obsdf.shape[0]))
    xe = np.empty((nobjects-1, obsdf.shape[0]))
    y = np.empty((nobjects-1, obsdf.shape[0]))
    ye = np.empty((nobjects-1, obsdf.shape[0]))
    
    colorcycle = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']
    #markercycle = ["x","+"]
    name_dict = runprops.get("names_dict")
    objectnames = []
    for i in name_dict.values():
        objectnames.append(i)
    #print(DeltaLong_Model)
    maxlong_index = np.argmax(DeltaLong_Model[1])
    maxlong = round(DeltaLong_Model[1][maxlong_index]+0.1,1)
    maxlat_index = np.argmax(DeltaLat_Model[1])
    maxlat = round(DeltaLat_Model[1][maxlat_index]+0.1,1)
    for j in range(20,len(DeltaLong_Model[0])-1):
        fig = plt.figure()
        plt.scatter(0,0)
        for i in range(1,nobjects):
            modelx[i-1,j] = DeltaLat_Model[i-1][j]
            modely[i-1,j] = DeltaLong_Model[i-1][j]
            
            x[i-1,:] = obsdf["DeltaLat_" + objectnames[i]].values
            xe[i-1,:] = obsdf["DeltaLat_" + objectnames[i] + "_err"].values
            y[i-1,:] = obsdf["DeltaLong_" + objectnames[i]].values
            ye[i-1,:] = obsdf["DeltaLong_" + objectnames[i] + "_err"].values
            
            plt.plot(DeltaLat_Model[i-1,(j-8):(j+1)], DeltaLong_Model[i-1,(j-8):(j+1)], colorcycle[i], alpha=0.3)
            plt.scatter(modelx[i-1,j], modely[i-1,j], color = colorcycle[i], label = objectnames[i],s=5)
            plt.errorbar(x[i-1,:], y[i-1,:], xerr = xe[i-1,:], yerr = ye[i-1,:], fmt = "ko", ms = 2, alpha=0.25)
    
        plt.axis('equal')
        plt.xlabel("Delta Latitude")
        plt.ylabel("Delta Longitude")
        plt.xlim(-maxlat,maxlat)
        plt.ylim(-maxlong,maxlong)
        plt.legend()
        fname = str(j).zfill(4)
        plt.savefig("spinny_astro/"+fname+"_astro.png")
        plt.close()
    
    print('fig ' + str(num) + ' done')


def auto_window(taus, c):
    m = np.arange(len(taus)) < c * taus
    if np.any(m):
        return np.argmin(m)
    return len(taus) - 1

def autocorr_new(y, c = 5.0):
    f = np.zeros(y.shape[1])
    for yy in y:
        f += emcee.autocorr.function_1d(yy)
    f /= len(y)
    taus = 2.0 * np.cumsum(f) - 1.0
    window = auto_window(taus, c)
    return taus[window]

def autocorrelation(sampler, objname, filename = "", thin = 1):
    # Getting chain for the first parameter to calculate different values
    chain = sampler.get_chain(thin = thin)
    
    nwalkers = sampler.nwalkers
    ndims = sampler.ndim
    nsteps = chain.shape[0]
    
    # Calculating values to calculate tau for
    # This chould be changed eventually
    N = np.exp(np.linspace(np.log(100), np.log(nsteps), 10)).astype(int)

    # Setting up array for tau estimates
    tau = np.empty( (len(N), ndims) )

    # Calculating tau for each value in N for each parameter
    for i in range(ndims):
        thischain = chain[:,:,i].T
        for j, n in enumerate(N):
            tau[j,i] = autocorr_new(thischain[:, :n])

    # Setting up to plot curves for tau in a grid
    x = 3
    y = ndims
    nrows = 1
    ncols = 3
    while x < y:
        y = y - x
        nrows += 1

    if ncols > ndims:
        ncols = ndims
    return np.mean(sampler.get_autocorr_time(quiet = True))

                
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

newpath = 'spinny_astro'
if not os.path.exists(newpath):
    os.makedirs(newpath)


fit_scale = pd.read_csv('fit_scale.csv',index_col=0)
float_names = runprops.get('float_names')
obsdf = pd.read_csv(objname+'_obs_df.csv',index_col=0)
geo_obj_pos = pd.read_csv('geocentric_'+objname+'_position.csv',index_col=0)
fixed_df = pd.read_csv('fixed_df.csv',index_col=0)
total_df_names = runprops.get('total_df_names')

#plots(backend, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names)
plot_astro(backend, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names)


'''   
    fig, ax = plt.subplots(nrows = nrows, ncols = ncols, sharey=True, 
                   gridspec_kw={'wspace': 0},
                   figsize = (6.4*(ncols),4.8*(nrows)), 
                   squeeze = False)
    fig.suptitle("Autocorrelation estimates")
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(False)
    plt.xlabel("number of samples, $N$")
    plt.ylabel(r"$\tau$ estimates")
    for i in range(nrows):
        for j in range(ncols):
            dim = i*ncols + j
            taus = ax[i,j].loglog(N, tau[:,dim], "o-", label="new")
            line = ax[i,j].plot(N, N / 50.0, "--k", label=r"$\tau = N/50$")
    fname = "../results/"+objname+"/autocorr" + filename + ".png"
    fig.savefig(fname, format='png')
'''
