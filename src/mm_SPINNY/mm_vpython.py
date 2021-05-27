import numpy as np
import pandas as pd
from mm_SPINNY.spinny_generate import *
from mm_SPINNY.quaternion import *
#from scipy.spatial.transform import Rotation as R
from vpython import *


def vpython_orbits(sys_df):
    
    t_arr = np.linspace(0,6*864000,2000)
    # default is integrate for 60 days. Change to compute length/step size based on object periods and velocities?
    T = len(t_arr)
    system = build_spinny(sys_df)

    N = system[0]
    names = system[1]
    phys_arr = system[2]
    orb_arr = system[3]
    spin_arr = system[4]
    quat_arr = system[5]
    
    spinny = evolve_spinny(N, names, phys_arr, orb_arr, spin_arr, quat_arr, t_arr)
    spinny_df = spinny[0]
    
    #create scene
    scene = canvas()
    
    temp_x = np.array([spinny_df["X_Pos_"+name].values.flatten() for name in names[1:]])
    temp_y = np.array([spinny_df["Y_Pos_"+name].values.flatten() for name in names[1:]])
    temp_z = np.array([spinny_df["Z_Pos_"+name].values.flatten() for name in names[1:]])
    
    i = 1
    for name in names[1:]:       
  
        globals()[name+'_x'] = prim2bary(phys_arr, temp_x)[i-1]
        globals()[name+'_y'] = prim2bary(phys_arr, temp_y)[i-1]
        globals()[name+'_z'] = prim2bary(phys_arr, temp_z)[i-1]
        
        globals()[name+'_obj'] = sphere(make_trail=True,trail_color = color.white, color=color.yellow) # create object
        
        globals()[name+'_obj'].radius = phys_arr[i,1]    
        i = i+1
        
    globals()[name+'_obj'].pos = vector(globals()[name+'_x'][0],globals()[name+'_y'][0],globals()[name+'_z'][0])
   
    j = 0
    scene.waitfor('keydown')
    while j < T:       
        rate(60)
        for name in names[1:]:
            globals()[name+'_obj'].pos = vector(globals()[name+'_x'][j],globals()[name+'_y'][j],globals()[name+'_z'][j])      
        j=j+1
        
def vpython_spins(sys_df):
    
    t_arr = np.linspace(0,6*864000,2000)
    # default is integrate for 60 days. Change to compute length/step size based on object periods and velocities?
    T = len(t_arr)
    system = build_spinny(sys_df)

    N = system[0]
    names = system[1]
    phys_arr = system[2]
    orb_arr = system[3]
    spin_arr = system[4]
    quat_arr = system[5]
    
    spinny = evolve_spinny(N, names, phys_arr, orb_arr, spin_arr, quat_arr, t_arr)
    spinny_df = spinny[0]
    
    #create scene
    scene = canvas()
    
    i = 1
    for name in names[1:]:
        globals()[name+'_inc'] = spinny_df["inclination_"+name].values.flatten()
        globals()[name+'_lan'] = spinny_df["longitude_ascending_"+name].values.flatten()
        globals()[name+'_aop'] = spinny_df["argument_periapse_"+name].values.flatten()
        globals()[name+'_mea'] = spinny_df["mean_anomaly_"+name].values.flatten()
       
        correct_lon = np.subtract(globals()[name+'_aop'],globals()[name+'_mea'])
        
        globals()[name+'_obl'] = np.subtract(spinny_df["obliquity_"+name].values.flatten(),globals()[name+'_inc'])
        globals()[name+'_prc'] = np.subtract(spinny_df["precession_"+name].values.flatten(),globals()[name+'_lan'])
        globals()[name+'_lon'] = np.subtract(spinny_df["longitude_"+name].values.flatten(),globals()[name+'_aop'])           
        
        globals()[name+'_obj'] = sphere()
        globals()[name+'_obj'].radius = phys_arr[i,1]
        
        # the primary object is built first at [0,0,0], the rest will be built based off its position/size)
        if i == 1:
            globals()[name+'_obj'].pos = vector(0.,0.,0.)
        else:
            globals()[name+'_obj'].pos = vector(globals()[names[i-1]+'_obj'].pos.x + (3.0 * globals()[names[i-1]+'_obj'].radius),0.,0.) 
        
        globals()[name+'_x'] = arrow(length = 2.0*globals()[name+'_obj'].radius, color = color.red)
        globals()[name+'_x'].pos = globals()[name+'_obj'].pos
        
        globals()[name+'_y'] = arrow(length = 2.0*globals()[name+'_obj'].radius, color = color.yellow)
        globals()[name+'_y'].pos = globals()[name+'_obj'].pos
        
        globals()[name+'_z'] = arrow(length = 2.0*globals()[name+'_obj'].radius, color = color.orange)
        globals()[name+'_z'].pos = globals()[name+'_obj'].pos
        i = i+1
    
    t = 0
    time_disp = label(text='Elapsed time:\n 0.00 days', pos = vector(0.,500.,0.), height = 10)
    while t <= T:
        rate(5)
        for name in names[1:]:
            phi =   globals()[name+'_prc'][t] 
            theta = globals()[name+'_obl'][t]
            psi =   globals()[name+'_lon'][t]
            #r = R.from_euler('zxz', [phi, theta, psi], degrees=True)
            
            rot_mat = quaternion(phi, theta, psi)[1] # convert the Euler angles to the associated rotation matrix
            
            xu = np.array([1.0, 0.0, 0.0])
            yu = np.array([0.0, 1.0, 0.0])
            zu = np.array([0.0, 0.0, 1.0])
            
            ax_x = xu.dot(rot_mat) #r.apply(xu)
            ax_y = yu.dot(rot_mat)
            ax_z = zu.dot(rot_mat) 
            
            #set/update the direction of each of the principle axes for each body
            globals()[name+'_x'].axis = 2.0*globals()[name+'_obj'].radius*vector(ax_x[0],ax_x[1],ax_x[2])
            globals()[name+'_y'].axis = 2.0*globals()[name+'_obj'].radius*vector(ax_y[0],ax_y[1],ax_y[2])
            globals()[name+'_z'].axis = 2.0*globals()[name+'_obj'].radius*vector(ax_z[0],ax_z[1],ax_z[2])
            time = str(t_arr[t]/(24.0*3600.0))
            time_disp.text = 'Elapsed time:\n'+time+' days'
            
        t = t+1    
        
def prim2bary(phys_arr, vec_arr):
    N = len(phys_arr)
    M = sum([phys_arr[n,0] for n in range(1,N)])
    T = len(vec_arr[0])
    
    bary_arr = np.empty(T)
    
    for t in range(0,T):
        bary_arr[t] = (M**(-1.0))*sum([phys_arr[n,0]*vec_arr[n-1,t] for n in range(1,N)])
        
    new_arr = np.array([np.subtract(vec_arr[n],bary_arr) for n in range(0,N-1)])
    
    return(new_arr)