
from mm_SPINNY.spinny import Spinny_System, Physical_Properties
from mm_SPINNY.spinny_plots import *
from mm_SPINNY.quaternion import *
import spiceypy as spice
import numpy as np
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from scipy.spatial.transform import Rotation as R

##### ALL UNITS SHOULD BE GIVEN IN km, sec, AND rad IN ORDER FOR THE COMPUTATIONS TO WORK

global G
G = 6.67e-20 # Gravitational constant in km
           
def orb2vec_ns(orb_arr,phys_arr,n): # converts orbital parameters to state vector 
    N = len(phys_arr)            # This MUST be done in an inertial frame or it cannot work

    a = orb_arr[0]
    e = orb_arr[1]
    w = orb_arr[2]
    i = orb_arr[3]
    O = orb_arr[4] 
    M = orb_arr[5]            
    
    list = [l for l in range(0,N)]
    mu = G*sum(phys_arr[list,0])
              
    if n == 0:   # scales the "orbit" of the primary according to the mass of the secondary
        orb = [a*(1-e),e,i,O,w,M,0.0,mu]
        cm = 1/(1 + (phys_arr[0,0]/phys_arr[1,0]))
        vec = -cm * spice.conics(orb,0)
        
    elif n == 1:   # scales the orbit or the secondary according to the mass of the primary
        orb = [a*(1-e),e,i,O,w,M,0.0,mu]
        cm = 1/(1 + (phys_arr[1,0]/phys_arr[0,0]))
        vec = cm * spice.conics(orb,0)    
        
    else:         # all the other objects just use the defined orbits
        orb = [a*(1-e),e,i,O,w,M,0.0,mu]
        vec = spice.conics(orb,0)

    return(vec) 


def vec2orb_ns(s,phys_objects,vec):  # converts a state vector to orbital parameters 
    N = len(phys_objects)
    
    mu = sum([phys_objects[n].mass for n in range(0,N)])
    orb = np.empty((T,8))
    for t in range(0,T):
        vec_t = vec[t,:]
        orb[t] = spice.oscelt(vec_t,s.t,mu)
            
    return(orb)               
                    # oscelt returns an array of 8 things:
                    # Perifocal distance.
                    # Eccentricity.
                    # Inclination.
                    # Longitude of the ascending node.
                    # Argument of periapsis.
                    # Mean anomaly at epoch.
                    # Epoch.
                    # Gravitational parameter.

                            
#### generate_system takes input from the dataframe and generates a 
#### Physical_Properties class and SPINNY object for each body in the system
    #Why is a plot command right here?
    #plot(s_df, names, runprops)
def generate_system_ns(N,name_arr,phys_arr,orb_arr,spin_arr,quat_arr,tolerance, runprops):
    
    # integration parameters
    #P = (2*np.pi)/np.sqrt(G*phys_arr[N-1,0]/(orb_arr[N-1,0]**3))

    tol = tolerance                         # integration tolerance
    h0P = 1.0                              # initial step size
    verbose = runprops.get("verbose")
    if verbose:
        print("Building SPINNY system...")
    s = Spinny_System(0.,h0=h0P,tol=tol)        # initializes object s, which is the SPINNY system
    
    for n in range(0,N):
        j2r20 = phys_arr[n,2]
        c22r20 = phys_arr[n,3]
        ax = phys_arr[n,1]

        # builds physical properties class for each object
        globals()['phys_'+str(n)] = Physical_Properties(G*phys_arr[n,0],"grav", J2R2=j2r20, C22R2=c22r20, c=ax)
        
        # creates an initial state vector for each body
        globals()['s_'+str(n)] = orb2vec_ns(orb_arr[n],phys_arr,n)           
    
        globals()['sp_rate_'+str(n)] = spin_arr[n,3]
        #r = R.from_euler('XZX', [spin_arr[n,0], spin_arr[n,1], spin_arr[n,2]])
        globals()['spin_'+str(n)] = np.array([0.0, 0.0, globals()['sp_rate_'+str(n)]])     
        
        s.add_object(name_arr[n],globals()['phys_'+str(n)],globals()['s_'+str(n)],globals()['spin_'+str(n)],quat_arr[n])
    
    s.move2bary()
    
    phys_objects = [globals()['phys_'+str(n)] for n in range(0,N)]
    
    return(s, phys_objects)
        
        
def spinny_integrate_ns(s, name_arr, phys_objects, t_arr, runprops): # evolves the SPINNY system to each time given in t_arr
    global T
    #epoch = 2453880.0 * 24.*3600.
    T = int(len(t_arr))
    N = int(len(s.arr0)/13)
    body_arr = np.empty((T,(N*6)))
    inertial_arr = np.empty((T,(N*13))) 
    spin_arr = np.empty((N,T,2))
    #spinvec isn't in spinny generate
    spinvec_arr = np.empty((N,T,3))
    quat_arr = np.empty((N,T,4))
    euler_arr = np.empty((N,T,3))
    hasspin_arr = np.ones(N, dtype=bool)
    L_arr = np.empty((N,T,3))
    E_arr = np.empty((N,T))
    verbose = runprops.get('verbose')
    if verbose:
        print("Evolving SPINNY...")

    for t in range(0,T):
        # Use s.get_state(n,0) with respect to the primary to ignore the motion of the Sun in our vectors,
        # but we take s.arr0 in order to get orbital parameters which might not make any sense in a primaricentric frame
        #print('t_arr', t_arr)
        s.evolve(t_arr[t]) #### <---- The actual integration
        #print("spinny_nosun Line 135")
        body_arr[t] = np.concatenate([s.get_state(n,0) for n in range(0,N)]) # taken with respect to the primary
        
        inertial_arr[t] = s.arr0
        
        for n in range(0,N):

            quat_n = s.get_quaternion(n)     # quaternion, (qr, qi, qj, qk)
            
            if quat_n.all() == 0.0:
                quat_n = [1.,0.,0.,0.]
                hasspin_arr[n] = False
            else:
                quat_n = quat_n
                hasspin_arr[n] = True
            
            qr = quat_n[0]
            qi = quat_n[1]
            qj = quat_n[2]
            qk = quat_n[3]
            
            quat_arr[n,t] = np.array([qi,qj,qk,qr])

            spinvec_arr[n,t,:] = s.get_spin(n)

    for t in range(0,T):
        for n in range(0,N):
            
            r_n = R.from_quat(quat_arr[n,t]) # convert quaternion to rotation
            obj_pole = r_n.apply([0.0,0.0,1.0],inverse=False) # get orientation of spin pole in world frame
            state = body_arr[t,(n*6):(n*6)+6] # get barycentric state vector for the body

            # These if statements should save some time in integration.
             # For the primary body, measure spin with respect to the second body's orbit so
             # that you don't get devide-by-zero warnings
            if n == 0 and hasspin_arr[n] == True:  
                state_sec = body_arr[t,((n+1)*6):((n+1)*6)+6]         
                h = np.cross(state_sec[:3],state_sec[3:])# Specific orbital angular momentum (of secondary)
                #print(h)
                orbit_pole = h/np.linalg.norm(h)  # compute direction of orbit normal
                spin_orbit_angle = np.arccos(np.dot(obj_pole,orbit_pole))*180./np.pi
                #print("here1")

            elif hasspin_arr[n] == False: # don't try to compute spin angles if spin isn't included
                spin_orbit_angle = 0.0
                #print("here2")

            else:
                h = np.cross(state[:3],state[3:]) # Specific orbital angular momentum
                orbit_pole = h/np.linalg.norm(h)  # compute direction of orbit normal
                #print(orbit_pole)
                spin_orbit_angle = np.arccos(np.dot(obj_pole,orbit_pole))*180./np.pi
                #print("here3")

            #print(spin_arr.shape)
            spin_arr[n,t,0] = spin_orbit_angle 
            spin_rate = np.linalg.norm(spinvec_arr[n,t,:]) # magnitude of spin vector
            spin_arr[n,t,1] = ((2.0*np.pi)/spin_rate) / 3600.0 # spin period in hours


            # converts quaternion to euler angles, using ZXZ rotation sequence   
            euler_arr[n,t] = r_n.as_euler('ZXZ') #quat2euler(quat_n) 


            # calculate angular momentum for the system to check for conservation
            I0 = phys_objects[n].I[0] # moment of inertia
            I1 = phys_objects[n].I[1] # moment of inertia
            I2 = phys_objects[n].I[2] # moment of inertia (this is the spin axis)
            w = spinvec_arr[n,t,:]

            L_body = np.array([I0 * w[0],I1 * w[1],I2 * w[2]]) 
            L_sp = r_n.apply(L_body,inverse=False) #translate angular momentum to world frame

            state_bary = s.get_state(n)
            h_bary = np.cross(state_bary[:3],state_bary[3:]) # specific orbital angular momentum, barycentric frame

            L_orb = phys_objects[n].mass * h_bary

            L_tot = np.add(L_orb, L_sp)

            L_arr[n,t] = L_tot

            # calculate total mechanical energy to check for conservation
            K_sp = np.sum(0.5 * np.array([I0 * w[0]**2.0,I1 * w[1]**2.0,I2 * w[2]**2.0]))
            K_orb = 0.5 * phys_objects[n].mass * np.linalg.norm(state_bary[3:])**2.0
            U_grav = np.sum([(phys_objects[n].mass * phys_objects[i].mass)/np.linalg.norm((s.get_state(n,i)[:3])) for i in range(0,N) if i is not n])
            E_tot = K_sp + K_orb - U_grav
            E_arr[n,t] = E_tot
    
    L_arr = np.sum(L_arr,axis=0)
    E_arr = np.sum(E_arr,axis=0)
    
    body_dict = {"Times":t_arr}
    if verbose:
        print("Constructing dataframe...")
    #print("spinny_no_sun line 222")
    for n in range(0,N):   
        body_dict.setdefault('X_Pos_'+name_arr[n] , body_arr[:,(n*6)+0])
        body_dict.setdefault('Y_Pos_'+name_arr[n] , body_arr[:,(n*6)+1])
        body_dict.setdefault('Z_Pos_'+name_arr[n] , body_arr[:,(n*6)+2])
        body_dict.setdefault('X_Vel_'+name_arr[n] , body_arr[:,(n*6)+3])
        body_dict.setdefault('Y_Vel_'+name_arr[n] , body_arr[:,(n*6)+4])
        body_dict.setdefault('Z_Vel_'+name_arr[n] , body_arr[:,(n*6)+5])
        
        cm1 = phys_objects[0].mass/(phys_objects[1].mass+phys_objects[0].mass)
        cm2 = phys_objects[1].mass/(phys_objects[1].mass+phys_objects[0].mass)
        
        vec = inertial_arr[:,(n*13):((n*13)+6)] - cm1*inertial_arr[:,0:6] - cm2*inertial_arr[:,13:13+6]
       
        ##### Attempt to subtract motion of the barycenter from state vectors
        #print([phys_objects[i].mass for i in range(1,N)])
        #vec_arr = np.array([body_arr[:,(i*6):((i*6)+6)] for i in range(0,N-1)])
        #M = sum([phys_objects[i].mass for i in range(1,N)])
        #T = len(vec_arr[0])
        #bary_arr = np.empty((T,6))
        
        #for t in range(0,T):
        #    bary_arr[t] = (M**(-1.0))*sum([phys_objects[i].mass*vec_arr[i-1,t] for i in range(1,N)])
        #    
        #new_vec = np.array([np.subtract(vec_arr[n,t],bary_arr[t]) for t in range(0,T)])
        
        orb_n = vec2orb_ns(s,phys_objects,vec)
            
        #semi-major axis calculated by a = (periapsis)/(1-e)
 
        body_dict.setdefault('semimajor_axis_'+name_arr[n]      , orb_n[:,0]/(1-orb_n[:,1]))
        body_dict.setdefault('eccentricity_'+name_arr[n]        , orb_n[:,1])
        body_dict.setdefault('inclination_'+name_arr[n]         , 180./np.pi*orb_n[:,2])
        body_dict.setdefault('longitude_ascending_'+name_arr[n] , 180./np.pi*orb_n[:,3])
        body_dict.setdefault('argument_periapse_'+name_arr[n]   , 180./np.pi*orb_n[:,4])
        body_dict.setdefault('mean_anomaly_'+name_arr[n]        , 180./np.pi*orb_n[:,5]) 
       
        # Euler angles, with respect to the orbit
        body_dict.setdefault('precession_'+name_arr[n], 180./np.pi*(euler_arr[n,:,0]) ) 
        body_dict.setdefault('obliquity_'+name_arr[n],  180./np.pi*(euler_arr[n,:,1]) )#-180./np.pi*orb_n[:,2] )
        body_dict.setdefault('longitude_'+name_arr[n],  180./np.pi*(euler_arr[n,:,2]) )
       
        body_dict.setdefault('spin_orbit_angle_'+name_arr[n], spin_arr[n,:,0])
        body_dict.setdefault('spin_period_'+name_arr[n],  spin_arr[n,:,1])
        
        #print(L_arr[:,0],L_arr[:,1],L_arr[:,2])
        body_dict.setdefault('Lx_'+name_arr[n], L_arr[:,0] )
        body_dict.setdefault('Ly_'+name_arr[n], L_arr[:,1] )
        body_dict.setdefault('Lz_'+name_arr[n], L_arr[:,2] )
        
        body_dict.setdefault('E_'+name_arr[n], E_arr )
                               
    spinny_df = pd.DataFrame(body_dict)
                               
    return(spinny_df)     

"""    
 useful things "s" contains:
    number of objects: nobj
    physical properties: phys
    names of bodies: called with s.idx("Name")
    the integration array: arr0 (x,y,z,vx,vy,vz,spinx,spiny,spinz,q1,q2,q3,q4) for EACH body
    do the objects have spin?: hasspin (boolean)
    time: t 
                     
"""                          

def build_spinny_ns(sys_df, runprops): 
    #print(sys_df)
    G = 6.67e-20 # Gravitational constant in km
    
    verbose = runprops.get('verbose')
    if verbose:
        print("Reading file to dataframe...")
    
    masses = -np.sort(-sys_df.loc["mass"].values.flatten()) # in DESCENDING ORDER of size, an array of the masses (sun is first)
    N = len(masses) 
    
    names = np.empty(N,dtype='object')
    #print(sys_df)
    for n in range(0,N):
        idx_n = int(np.where(sys_df==masses[n])[1].flatten())
        name_n = (sys_df.columns[idx_n])
        names[n] = name_n
    
    phys_arr = np.zeros((N,4))
    orb_arr = np.zeros((N,6))
    spin_arr = np.zeros((N,4))
    quat_arr = np.zeros((N,4))
    
    i = 0
    for name in names: # for each body in the system, added in order of descending mass:
        
        if sys_df.loc["mass",name] != 0.00:
            mass_n = sys_df.loc["mass",name] 
            
        else:
            mass_n = 0.0 
            
        if sys_df.loc["axis",name] != 0.00:
            ax_n = sys_df.loc["axis",name]
        else:
            ax_n = 1.0
            
        if sys_df.loc["j2r2",name] != 0.00:
            j2r2_n = sys_df.loc["j2r2",name]
        else:
            j2r2_n = 0.0    
        
        if sys_df.loc["c22r2",name] != 0.00:
            c22r2_n = sys_df.loc["c22r2",name]
        else:
            c22r2_n = 0.0    
                
        # set physical properties array
        phys_arr[i] = np.array([mass_n, ax_n, j2r2_n, c22r2_n])
        
        if sys_df.loc["sma",name] != 0.00:
            sma_n = sys_df.loc["sma",name]
        else:
            sma_n = 0.0
        
        if sys_df.loc["ecc",name] != 0.00:
            ecc_n = sys_df.loc["ecc",name]
        else:
            ecc_n = 0.0
            
        if sys_df.loc["aop",name] != 0.00:
            aop_n = sys_df.loc["aop",name]
        else:
            aop_n = 0.0
            
        if sys_df.loc["inc",name] != 0.00:
            inc_n = sys_df.loc["inc",name]
        else:
            inc_n = 0.0 
            
        if sys_df.loc["lan",name] != 0.00:
            lan_n = sys_df.loc["lan",name]
        else:
            lan_n = 0.0
        
        if sys_df.loc["mea",name] != 0.00:
            mea_n = sys_df.loc["mea",name]
        else:
            mea_n = 0.0
            
        # set orbital properties array
        orb_arr[i] = np.array([sma_n, ecc_n, aop_n, inc_n, lan_n, mea_n]) 
        i = i + 1
        
    # set orbit of primary equal to orbit of secondary (it will be scaled later)

    orb_arr[0] = orb_arr[1]
    i = 0
    
    for name in names:
        # precession, obliquity angles are measured with respect to the ECLIPTIC, not the body's orbit
        # default values are set to be aligned with the orbit (LAN for prec, inc for obliq, AOP for longitude)
        
        if sys_df.loc["sp_prc",name] != 0.00:
            sp_prc_n = sys_df.loc["sp_prc",name]
        else:
            sp_prc_n = orb_arr[i,4]
            
        if sys_df.loc["sp_obl",name] != 0.00:
            sp_obl_n = sys_df.loc["sp_obl",name]
        else:
            sp_obl_n = orb_arr[i,3]
            
        if sys_df.loc["sp_lon",name] != 0.00:
            sp_lon_n = sys_df.loc["sp_lon",name]
        else:
            sp_lon_n = orb_arr[i,2]
            
        if sys_df.loc["sp_rate",name] != 0.00:
            sp_rate_n = sys_df.loc["sp_rate",name]
        else:
            # the default is equal to the orbital period, calculated from semi-major axis, sum of masses of body and primary 
            sp_rate_n = (orb_arr[i,0]**3/(G*(phys_arr[i,0]+phys_arr[0,0])) )**(-0.5) 
            
        # set spin properties array
        spin_arr[i] = np.array([sp_prc_n, sp_obl_n, sp_lon_n, sp_rate_n])

        # create orientation quaternions for each body from spin/orbit data
        # generate quaternion for obliquity angle of rotation axis
        # Euler angles: precession, obliquity, longitude, relative to the ecliptic
        
        
        r = R.from_euler('ZXZ', [spin_arr[i,0], spin_arr[i,1], spin_arr[i,2]])
        quat_i = r.as_quat() 
        
        qi = quat_i[0]
        qj = quat_i[1]
        qk = quat_i[2]
        qr = quat_i[3]
        # set quaternion array
        # the order of the quaternion must be changed, since SPINNY needs scalar-first form
        quat_arr[i] = np.array([qr,qi,qj,qk]) #obliq_quat_n[0]
        
        i = i+1

    return(N, names, phys_arr, orb_arr, spin_arr, quat_arr)
    
    
def evolve_spinny_ns(N, names, phys_arr, orb_arr, spin_arr, quat_arr, t_arr, tol, runprops):
    
    start_time = time.time()
    s = generate_system_ns(N,names,phys_arr,orb_arr,spin_arr,quat_arr,tol, runprops)
    
    spinny = s[0]
    phys_arr = s[1]
    
    s_df = spinny_integrate_ns(spinny, names, phys_arr,t_arr, runprops)
    verbose = runprops.get('verbose')
    if verbose:
        print("Done.")

    seconds_time = time.time() - start_time

    if seconds_time >= 60.0:
        minutes_time = seconds_time/60.0
        if verbose:
            print("Completed in "+str(minutes_time)+" minutes.")

    else:
        if verbose:
            print("Completed in "+str(seconds_time)+" seconds.")

    if verbose:
        print("")
    np.insert(names,0,"Sun")
    return(s_df, phys_arr, names)


    
    
    
