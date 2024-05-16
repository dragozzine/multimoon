import numpy as np
import pandas as pd
import sys
from astropy import units as u
from astropy.coordinates import SkyCoord

from astropy.coordinates import GeocentricTrueEcliptic
sys.path.insert(1, 'mm_SPINNY/')
sys.path.insert(2, 'prep/')
from mm_SPINNY.spinny_vector import generate_vector
import prep.latlon_transform as latlon_transform
import mm_priors as prior
import mm_param
import random
import commentjson as json

from astropy.time import Time
import mm_relast
from csv import writer
import latlon_transform
import os
import time
from scipy.stats import chi2
from astroquery.jplhorizons import Horizons
import mpmath as mp
#from func_timeout import func_timeout, FunctionTimedOut


class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

"""
Inputs:
1)The Parameters dataframe
2) The Observation Dataframe
Outputs:
1) The chi-squared number of the likelihood
"""
# calculates the chi-square for parameters given observations
def mm_chisquare(paramdf, obsdf, runprops, geo_obj_pos, gensynth = False):

    #print('paramdf',paramdf)
    numObj = runprops.get("numobjects")
    verbose = runprops.get("verbose")
    pd.set_option('display.max_columns', None)
    names = []
    for i in range(1,numObj+1):
        names.append('name_'+str(i))
        if not 'name_'+str(i) in paramdf.columns:
            print('The parameter name_' + str(i)+ ' is not found in the parameter dataframe.')
            sys.exit()
        
    
    # use parameters dataframe with one set of parameters and observation times to call SPINNY to get the model dataframe
    # SPINNY returns (full?) state vector of each object (in primaricentric ecliptic J2000 coordinates) at each time
    
    # parameters dataframe "paramdf" is 1 row and many columns
    # N objects, each of which is a point-mass, oblate, or triaxial
       # point-masses: mass, 6 orbital elements
       # oblate: + J2R2, spinc, splan
       # triaxial: + C22R2, spaop, sprate
       # dynamicstoincludeflags='2010'
    # Include Sun or not
    # it has already been "checked" to make sure that all of the columns that are needed are there
    # mass_0, sma_0, ecc_0, ..., mea_0 are the orbital elements of the Sun 
       # at the epoch in primaricentric coordinates
    # when columns are not needed they are typically absent

    # MultiMoon units are km, kg, deg, seconds

    obsdf = obsdf.sort_values(by=['time'])
    #time_arr = np.sort(obsdf['time'].values.flatten())
    #print('OBSDF: ',obsdf)
    time_arr = obsdf['time'].values.flatten()# gets an array of observation times from the obs dataframe

    # Setting times relative to the epoch
    epoch = runprops.get("epoch_SJD")
    time_arr = time_arr - epoch

    # Sorts them into ascending order
#    import logging
    #print('time_Arr', time_arr)
    #print('obsdf', obsdf)

    begin = time.time()
    #print(paramdf)
    #print('param_df:', paramdf)
    try:
        time_arr_sec = time_arr*86400
        #vec_df = func_timeout(5,generate_vector,args=(paramdf, time_arr_sec, runprops))
        #print(paramdf)
        #print('time_arr_sec: ',time_arr_sec)
        vec_df = generate_vector(paramdf, time_arr_sec, runprops)
        
    #except FunctionTimedOut:
    #    print('Spinny took longer than 5 seconds to run 1 walker-step:\n')
    #    return np.inf
#    except Exception as e:
#        logging.exception('')
#        return np.inf
    except Exception as e:
        print('There was an error thrown within spinny:\n', e)
        rows = obsdf.shape[0]
        return -np.inf, np.ones(((numObj-1)*2, rows))*10000
    names_dict = runprops.get("names_dict")
    names=[0 for i in range(numObj)]
    for i in range(0,numObj):
        names[i] = names_dict.get("name_"+str(i+1))
    end = time.time()
    #print(end-begin,' seconds')
    #if end-begin > 2:
    #    print('A step took ', end-begin,' seconds')
    #    print(paramdf)

    # vec_df is a dataframe with len(time_arr) rows and
    # columns are state parameters x nobjects
    # Example: vecdf["X_Pos_"+paramsdf["name_2"]] gets the x position of object 2
    # ecliptic (J2000) coordinates
    # km, kg, rad, s
    # primaricentric 

    name_1 = "X_Pos_"+names[0]
    #print(vec_df)
    if (vec_df[name_1][0] != 0.0):
        print("Not primaricentric like I thought!")
        rows = obsdf.shape[0]
        return -np.inf, np.ones(((numObj-1)*2, rows))*10000
        #print("vec_df[name_1] = ", vec_df)
    
    Model_DeltaLong = np.zeros((numObj-1,len(time_arr)))
    Model_DeltaLat = np.zeros((numObj-1,len(time_arr)))
    if runprops.get('includesun') == 1:
        #print(vec_df)
        vec_df = vec_df.drop(['X_Pos_Sun', 'Y_Pos_Sun', 'Z_Pos_Sun', 'X_Vel_Sun', 'Y_Vel_Sun', 'Z_Vel_Sun'], axis=1)

    positionData = np.zeros((numObj*3,len(time_arr)))
        
    for i in range(0,numObj):
        positionData[3*i] = vec_df["X_Pos_"+names[i]]
        positionData[3*i+1] = vec_df["Y_Pos_"+names[i]]
        positionData[3*i+2] = vec_df["Z_Pos_"+names[i]]

        # tind = index/row number of vec_df corresponding to this time
        
        # for object from 2:N
             # gather relative positions
             # thisdx = vec_df[tind,"X_Pos_"+paramdf["name_"+str(object)] [ - X_Pos_1=0]
             # same for dy and dz

        # Implicitly assume that observer is close to geocenter (within ~0.01 AU)
        
        # obs_to_prim_pos = vectors of observer to primary
        # prim_to_sat__pos = vectors of primary to satellite

    obs_to_prim_pos = [positionData[0]+geo_obj_pos['x'].tolist(),positionData[1]+geo_obj_pos['y'].tolist(),positionData[2]+geo_obj_pos['z'].tolist()]
    
    prim_to_sat_pos = [positionData[1*3],positionData[1*3+1],positionData[1*3+2]]
    Model_DeltaLong, Model_DeltaLat = mm_relast.convert_ecl_rel_pos_to_geo_rel_ast(obs_to_prim_pos, prim_to_sat_pos)
    return obs_to_prim_pos,Model_DeltaLong,Model_DeltaLat


def pos_draw(sampler, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names):
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
    #print(sampler)    
    
    for i in range(runprops.get('numobjects')-1):
        #print(i,float_names)        
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
    if isinstance(runprops.get('thin_plots'), int):          
        thin_plots = runprops.get('thin_plots')
        if runprops.get('thin_run') and thin_plots > 1:
            print('Warning: You thinned your chain as you ran the data, and have a thin_plots parameter > 1, this could lead to an empty chain being read in.')
        #burnin = burnin/thin_plots
        #clusterburn = clusterburn/thin_plots     
    else:
        thin_plots = 1
        
        
#    chain = sampler.get_chain(discard=int(burnin+clusterburn),flat = False, thin=thin_plots) 
    #print(os.getcwd(), '')
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
    shortchain = sampler.get_chain(flat = False, thin= thin_plots)
    numparams = shortchain.shape[2]
    numwalkers = shortchain.shape[1]
    numgens = shortchain.shape[0]
    del shortchain
 
    # Take chain "fit" values and make them into real values
    for i in range(numparams):
        chain[:,:,i] = chain[:,:,i]*fit[i]
    #print('First chain: ',chain.shape)

    
    fitparam_chain = np.zeros((1,numwalkers,numgens))

    #print(fitparam_chain.shape)    
    fitparam_names = []    
    # Now de-transform the chain
    #print("Starting un transformations")
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
                #'''
                inc_new = chain[:,:,int(inc_lan_index[b*2])]
                lan_new = chain[:,:,int(inc_lan_index[b*2+1])]
                fitparam_chain = np.concatenate((fitparam_chain, np.array([inc_new.T])),axis=0)
                fitparam_chain = np.concatenate((fitparam_chain, np.array([lan_new.T])),axis=0)
                fitparam_names.append('equinoctial_q_'+str(b+1))
                fitparam_names.append('equinoctial_p_'+str(b+1))
                lan = (np.arctan2(inc_new,lan_new)*180/np.pi)%360
                chain[:,:,int(inc_lan_index[b*2+1])] = lan
                inc = (np.arccos(inc_new/np.sin(lan*np.pi/180))*2*180/np.pi)%180
                chain[:,:,int(inc_lan_index[b*2])] = inc                
                
                
                '''
                inc_new = chain[:,:,int(inc_lan_index[b*2])]
                lan_new = chain[:,:,int(inc_lan_index[b*2+1])]
                fitparam_chain = np.concatenate((fitparam_chain, np.array([inc_new.T])),axis=0)
                fitparam_chain = np.concatenate((fitparam_chain, np.array([lan_new.T])),axis=0)
                fitparam_names.append('equinoctial_q_'+str(b+1))
                fitparam_names.append('equinoctial_p_'+str(b+1))
                lan = (np.arctan2(inc_new,lan_new)*180/np.pi)%360
                chain[:,:,int(inc_lan_index[b*2+1])] = lan
                inc = (np.arctan2(inc_new,np.sin(lan*np.pi/180))*2*180/np.pi)%180
                chain[:,:,int(inc_lan_index[b*2])] = inc#'''
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
                spinc = (np.arccos(spinc_new/np.sin(splan*np.pi/180))*2*180/np.pi)%180
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
    fitparam_chain = fitparam_chain[int(burnin/thin_plots+clusterburn/thin_plots) :: 1]
        
    #print("Un transforming done")

    # Cutting up chain
    full_chain = np.copy(chain)
    chain = chain[int(burnin/thin_plots+clusterburn/thin_plots) :: 1]
    rand1 = np.random.random_integers(0,len(chain))
    rand2 = np.random.random_integers(0,len(chain[0]))
    #print(chain[rand1,rand2])
    #print(float_names)    
    return chain[rand1,rand2]


def run_eph_pos(k,paramdf, backend, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names):
    pdraw = pos_draw(backend, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names)
    
    epoch = runprops.get('epoch_SJD')
    time_to_L2 = 5.016722408026756/24/60/60
    
    
    #best_llhood = pd.read_csv(filename+'/best_likelihoods.csv')
    #best=best_llhood.iloc[-1]

    paramdf = paramdf.T
    paramdf['name_1'] = objname
    paramdf['name_2'] = runprops.get('names_dict').get('name_2')
    paramdf['time'] = epoch
    paramdf['mass_1']['mean'] = pdraw[0]
    paramdf['mass_2']['mean'] = pdraw[1]
    paramdf['sma_2']['mean'] = pdraw[2]
    paramdf['ecc_2']['mean'] = pdraw[3]
    paramdf['inc_2']['mean'] = pdraw[5]
    paramdf['aop_2']['mean'] = pdraw[4]
    paramdf['lan_2']['mean'] = pdraw[6]
    paramdf['mea_2']['mean'] = pdraw[7]
    
    total_model = pd.DataFrame()
    total_ra = np.zeros(150*365)
    total_dec = np.zeros(150*365)
    total_dist = np.zeros(150*365)
    #print(paramdf)
    for i in range(150):
        #if i%10 == 0:
            #print(i)
        times = np.arange(2455680+i*36.5,2455680+36.5*(i+1),0.1)
        ourKBO = Horizons(id=objname,location=location,epochs = times)
        ephKBO = ourKBO.ephemerides()['RA','DEC','datetime_jd']
        vecKBO = ourKBO.vectors(aberrations = 'astrometric')['lighttime','x','y','z']
    
        jdTime= ephKBO['datetime_jd']
        lightTime = vecKBO['lighttime']
        kboTime=jdTime-lightTime
    
        geo_obj_pos = pd.DataFrame({'kboTIME':kboTime,'x':vecKBO['x']*149597870.7 ,'y':vecKBO['y']*149597870.7 ,'z':vecKBO['z']*149597870.7 })
        geo_obj_pos['time'] = jdTime-time_to_L2
    
        data = np.array(mm_chisquare(paramdf, geo_obj_pos, runprops, geo_obj_pos))
        total_ra[365*i:365*(i+1)] = np.array(ephKBO['RA'])
        total_dec[365*i:365*(i+1)] = np.array(ephKBO['DEC'])
    
    #print(data[0])
    
        #print(paramdf)
    
        radec_data = np.array([ephKBO['RA'],ephKBO['DEC']])
    
        dist = ourKBO.vectors(aberrations = 'astrometric')['range']
    
        total_dist[365*i:365*(i+1)] = np.array(dist)
        
        dateList = []
        for j in jdTime:
            jd = Time(j,format='jd')
            dateList.append(jd)
    
        primC = SkyCoord(ra=ephKBO['RA'], dec=ephKBO['DEC'], frame='gcrs', obstime = dateList, distance = dist)
        primEcl = primC.transform_to(GeocentricTrueEcliptic(equinox='J2000'))
        
        Lat_Prim = primEcl.lat.degree
        Long_Prim = primEcl.lon.degree
    
    #latlon = latlon_transform.convert_to_primary_centric(paramdf, ['Eris','Dysnomia'], 2, 500)
    #print(latlon)
    
        vector = pd.DataFrame(np.array(data[0]).T,columns=['x','y','z'])
    #print(vector)
    #print(data[1],data[2])
        Model_df = pd.DataFrame(np.array([geo_obj_pos['time'],Lat_Prim,Long_Prim,data[1],data[2]]).T,columns=['time','Lat_prim','Lon_prim','DeltaLat','DeltaLon'])
        total_model = pd.concat([total_model, Model_df])
        
        
    from astropy.coordinates import GCRS
    
    arr = np.arange(0,365*150)
    ind = pd.Index(arr)
    
    total_model.index = ind
    from astropy.coordinates import GCRS
    dateList = total_model['time']
    
    for i in range(len(dateList)):
        jd = Time(dateList[i],format='jd')
        dateList[i]=jd
    total_model['Delta-RA_Secondary'] = np.zeros(len(total_model))
    total_model['Delta-DEC_Secondary'] = np.zeros(len(total_model))
    total_model['Delta-RA_Secondary-err'] = np.zeros(len(total_model))
    total_model['Delta-DEC_Secondary-err'] = np.zeros(len(total_model))
    
    total_model['RA_prim'] = total_ra
    total_model['DEC_prim'] = total_dec
    #time,Delta-RA_Secondary,Delta-RA_Secondary-err,Delta-DEC_Secondary,Delta-DEC_Secondary-err
    for j in range(len(total_model)):
        coord_sky = SkyCoord((total_model['DeltaLat'][j]/3600+total_model['Lat_prim'][j])*u.degree, (total_model['DeltaLon'][j]/3600+total_model['Lon_prim'][j])*u.degree, frame='geocentrictrueecliptic', obstime = dateList[j], distance = total_dist[j]*u.AU,unit=(u.deg,u.deg))
        prim_co = SkyCoord((total_model['Lat_prim'][j])*u.degree, (total_model['Lon_prim'][j])*u.degree, frame='geocentrictrueecliptic', obstime = dateList[j], distance = total_dist[j]*u.AU,unit=(u.deg,u.deg))
        
        coord_sky = coord_sky.transform_to(GCRS())
        prim_co = prim_co.transform_to(GCRS())
    
        total_model['Delta-RA_Secondary'][j] = -(prim_co.ra.degree-coord_sky.ra.degree)*3600
        total_model['Delta-DEC_Secondary'][j] = -(prim_co.dec.degree-coord_sky.dec.degree)*3600
        
    total_model.to_csv('../results/Eris/David_Ephemerides/'+str(k)+'_posdraw.csv')

import emcee
objname = 'Eris'
file = 'Eris_2023-09-14_21.42.55_00_robust_2018_v6'
filename = '../results/'+objname+'/'+file
getData = ReadJson(filename+'/runprops.txt')
runprops = getData.outProps()
obsdf = pd.read_csv(filename+'/'+objname+'_obs_df.csv')
location = '@SEMB-L2'
backend = emcee.backends.HDFBackend(filename+'/chain.h5')
getData = ReadJson(filename+'/runprops.txt')

runprops = getData.outProps()

fit_scale = pd.read_csv(filename+'/fit_scale.csv',index_col=0)
float_names = runprops.get('float_names')
obsdf = pd.read_csv(filename+'/'+objname+'_obs_df.csv',index_col=0)
geo_obj_pos = pd.read_csv(filename+'/geocentric_'+objname+'_position.csv',index_col=0)
fixed_df = pd.read_csv(filename+'/fixed_df.csv',index_col=0)
total_df_names = runprops.get('total_df_names')

paramdf = pd.read_csv(filename+'/'+objname+'_init_guess.csv',index_col=0)

for k in range(100):
    run_eph_pos(k,paramdf, backend, fit_scale, float_names, obsdf, runprops, geo_obj_pos, fixed_df, total_df_names)
    
    