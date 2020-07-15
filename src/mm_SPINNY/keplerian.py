import numpy as np
import spiceypy as spice
import pandas as pd
from spinny_generate import * 
from time import ctime
import time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import sys
sys.path.append("..")
import mm_runprops
runprops = mm_runprops.runprops

# computes a two-body, Keplerian integration using spiceypy
def kepler_2body(sys_df,t_arr):
    
    sys_df = sys_df.iloc[0]
    
    runprops = mm_runprops.runprops
    verbose = runprops.get("verbose")
    
    G = 6.674e-20 # Gravitational constant in km
    T = len(t_arr)
    
    names = [sys_df["name_1"], sys_df["name_2"]]  
    masses = sys_df["mass_1"] + sys_df["mass_2"]   
    
    mu = G*masses # compute the gravitational parameter for the system
    
    a2 = sys_df["sma_2"]
    e2 = sys_df["ecc_2"]
    i2 = sys_df["inc_2"]*(np.pi/180.0)
    O2 = sys_df["lan_2"]*(np.pi/180.0)
    w2 = sys_df["aop_2"]*(np.pi/180.0)
    M2 = sys_df["mea_2"]*(np.pi/180.0)
    p2 = a2*(1-e2) # perifocal distance
    
    T0 = t_arr[0] # epoch (seconds past J2000)
    
    orb2 = [p2,e2,i2,O2,w2,M2,T0,mu]
    
    vec_arr = np.empty((T,6))
    orb_arr = np.empty((T,8))
  
    start_time = time.time()
    
    if verbose:
        print("Running two-body Keplerian integration...")
    
    for t in range(0,T):
    
        vec_arr[t] = spice.conics(orb2,t_arr[t])
        orb_arr[t] = spice.oscelt(vec_arr[t],t_arr[t],mu)
                         
        # oscelt returns an array of 8 things:
        # Perifocal distance.
        # Eccentricity.
        # Inclination.
        # Longitude of the ascending node.
        # Argument of periapsis.
        # Mean anomaly at epoch.
        # Epoch.
        # Gravitational parameter.
    if verbose:
        print("Done.")

        seconds_time = time.time() - start_time

        if seconds_time >= 60.0:
            minutes_time = seconds_time/60.0
            print("Completed in "+str(minutes_time)+" minutes.")

        else:
            print("Completed in "+str(seconds_time)+" seconds.")

    body_dict = {"Times":t_arr}
    if verbose:
        print("Constructing dataframe...")
    
   
    body_dict.setdefault('X_Pos_'+names[1] , vec_arr[:,0])
    body_dict.setdefault('Y_Pos_'+names[1] , vec_arr[:,1])
    body_dict.setdefault('Z_Pos_'+names[1] , vec_arr[:,2])
    body_dict.setdefault('X_Vel_'+names[1] , vec_arr[:,3])
    body_dict.setdefault('Y_Vel_'+names[1] , vec_arr[:,4])
    body_dict.setdefault('Z_Vel_'+names[1] , vec_arr[:,5])

    body_dict.setdefault('sma_'+names[1] , orb_arr[:,0]/(1-orb_arr[:,1]))
    body_dict.setdefault('ecc_'+names[1] , orb_arr[:,1])
    body_dict.setdefault('inc_'+names[1] , orb_arr[:,2]*(180.0/np.pi))
    body_dict.setdefault('lan_'+names[1] , orb_arr[:,3]*(180.0/np.pi))
    body_dict.setdefault('aop_'+names[1] , orb_arr[:,4]*(180.0/np.pi))
    body_dict.setdefault('mea_'+names[1] , orb_arr[:,5]*(180.0/np.pi))
    
    body_dict.setdefault('X_Pos_'+names[0] , 0.0)
    body_dict.setdefault('Y_Pos_'+names[0] , 0.0)
    body_dict.setdefault('Z_Pos_'+names[0] , 0.0)
    body_dict.setdefault('X_Vel_'+names[0] , 0.0)
    body_dict.setdefault('Y_Vel_'+names[0] , 0.0)
    body_dict.setdefault('Z_Vel_'+names[0] , 0.0)
    
    kepler_df = pd.DataFrame(body_dict)
    
    return(kepler_df, names)
    
def kepler_nbody(sys_df,t_arr): # runs Keplerian integrations for systems with ONE massive body and N massless bodies
        
    sys_df = sys_df.iloc[0]
    
    runprops = mm_runprops.runprops
    verbose = runprops.get("verbose")
    
    G = 6.674e-20 # Gravitational constant in km
    T = len(t_arr)
    
    N = runprops.get("numobjects")

    names = [sys_df["name_"+str(n)] for n in range(1,N)]
    mass_1 = sys_df["mass_1"]

    mu = G*mass_1 # compute the gravitational parameter for the system (only one body)
    
    T0 = t_arr[0] # epoch (seconds past J2000)

    orb_init = np.empty((N-1,8))
    vec_init = np.empty((N-1,6))
    for n in range(2,N): # exclude the primary, it should not have any orbital elements given
        a_n = sys_df["sma"+str(n)].iloc[0]
        e_n = sys_df["ecc"+str(n)].iloc[0]
        i_n = sys_df["inc"+str(n)].iloc[0]
        O_n = sys_df["lan"+str(n)].iloc[0]
        w_n = sys_df["aop"+str(n)].iloc[0]
        M_n = sys_df["mea"+str(n)].iloc[0]
        p_n = a_n*(1-e_n) # periapsis distance
    
        orb_init[n] = [p_n,e_n,i_n,O_n,w_n,M_n,T0,mu] # set orbital params list
     
        vec_init[n] = spice.conics(orb_init[n],t_arr[0])
    
        start_time = time.time()
    
    orb_arr = np.empty((N-1,T,8))
    vec_arr = np.empty((N-1,T,6))
    print("Running "+str(N)+"-body Keplerian integration...")
    for t in range(0,T):
        for n in range(1,N):
            vec_arr[n,t] = spice.conics(orb_init[n],t_arr[t])
            orb_arr[n,t] = spice.oscelt(vec_arr[t],t,mu)
                         
            # oscelt returns an array of 8 things:
            # Perifocal distance.
            # Eccentricity.
            # Inclination.
            # Longitude of the ascending node.
            # Argument of periapsis.
            # Mean anomaly at epoch.
            # Epoch.
            # Gravitational parameter.
    print("Done.")

    seconds_time = time.time() - start_time

    if seconds_time >= 60.0:
        minutes_time = seconds_time/60.0
        print("Completed in "+str(minutes_time)+" minutes.")

    else:
        print("Completed in "+str(seconds_time)+" seconds.")
        
    body_dict = {"Times":t_arr}
    print("Constructing dataframe...")
    
    for name in names[1:]:
        body_dict.setdefault('X_Pos_'+name , vec_arr[:,0])
        body_dict.setdefault('Y_Pos_'+name , vec_arr[:,1])
        body_dict.setdefault('Z_Pos_'+name , vec_arr[:,2])
        body_dict.setdefault('X_Vel_'+name , vec_arr[:,3])
        body_dict.setdefault('Y_Vel_'+name , vec_arr[:,4])
        body_dict.setdefault('Z_Vel_'+name , vec_arr[:,5])

        body_dict.setdefault('sma_'+name , orb_arr[:,0]/(1-orb_arr[:,1]))
        body_dict.setdefault('ecc_'+name , orb_arr[:,1])
        body_dict.setdefault('inc_'+name , orb_arr[:,2])
        body_dict.setdefault('lan_'+name , orb_arr[:,3])
        body_dict.setdefault('aop_'+name , orb_arr[:,4])
        body_dict.setdefault('mea_'+name , orb_arr[:,5])

    kepler_df = pd.DataFrame(body_dict)
    
    return(kepler_df, names)
    
def kepler_save(kepler_df, names):   
    save_yn = str(raw_input("Do you want to save this data? (Y/N): "))
    save_yn = save_yn.upper()
    
    if save_yn=="Y":
        print("Generating .csv...")
        t_current = ctime().replace(" ","_")
        filename = names[0]+"_KEPLER_"+t_current+".csv"
        kepler_df.to_csv(filename)
        print("Data from SpiceyPy saved to the local file as "+filename)
        kepler_plot_q(kepler_df,names)

    elif save_yn == "N":
        print("")
        kepler_plot_q(kepler_df,names)
    else:
        print("")
        print('Invalid Response.')
        return save(kepler_df,names) 
    
    
    
def kepler_plot_q(kepler_df,names): 
    plot_yn = str(raw_input("Do you want me to generate figures from these data? (Y/N): "))
    plot_yn = plot_yn.upper()
    
    if plot_yn == "Y":
        kepler_plot(kepler_df,names)
    elif plot_yn == "N":
        print("\n Returning to main menu...")
        return main_menu()
    else:
        print("")
        print('Invalid Response.')
        return plot_q(kepler_df,names)    

def kepler_plot(plot_df,names):
   
    t_arr = plot_df['Times'].values.flatten()
    
    globals()[names[1]+'_x'] = plot_df["X_Pos_"+names[1]]
    globals()[names[1]+'_y'] = plot_df["Y_Pos_"+names[1]]
    globals()[names[1]+'_z'] = plot_df["Z_Pos_"+names[1]]
    
    fig2,ax = plt.subplots(figsize=(10,10))
    ax2 = plt.subplot2grid((4, 4), [0, 1], 2 , 2)           
    ax3 = plt.subplot2grid((4, 4), (2, 2), 2 , 2, sharey=ax2)
    ax4 = plt.subplot2grid((4, 4), (2, 0), 2 , 2, sharey=ax2)
   
    ax2.set_aspect('equal')
    ax3.set_aspect('equal')
    ax4.set_aspect('equal')
    
    ax2.grid()
    ax3.grid()
    ax4.grid()

    ax2.set_xlabel('x (kilometers)', fontsize = 18)
    ax2.set_ylabel('z (kilometers)', fontsize = 18)

    ax3.set_xlabel('z (kilometers)', fontsize = 18)
    ax3.set_ylabel('y (kilometers)', fontsize = 18)
    plt.setp(ax3.get_yticklabels(), visible=False)
    
    ax4.set_xlabel('x (kilometers)', fontsize = 18)
    ax4.set_ylabel('y (kilometers)', fontsize = 18)
    
    t_current = ctime().replace(" ","_")
    filename = names[0]+"_KEPLER_figures_"+t_current+".pdf"
    
    with PdfPages(filename) as pdf:
        t = t_arr/(3600*24)        
        if t[-1] > 1000:
            t = t/365.25
            int_time = str(t[-1])+" Years" 
        else:
            int_time = '%.3f Days' % t[-1]
            
        ax2.plot([0.0],[0.0], marker='o', color='black',label=names[0])
        ax3.plot([0.0],[0.0], marker='o', color='black',label=names[0])
        ax4.plot([0.0],[0.0], marker='o', color='black',label=names[0])
        fig2.suptitle("Primaricentric Orbit of "+names[1]+" for "+int_time ,fontsize=18, fontweight='bold')
        
        nx = globals()[names[1]+'_x']
        ny = globals()[names[1]+'_y']
        nz = globals()[names[1]+'_z']

        if max(nx) > 10000.0: # adjust the units of orbit plots, for readability
            nx = nx/1000.0
            ny = ny/1000.0
            nz = nz/1000.0

            ax2.set_xlabel('x ($10^3$ kilometers)', fontsize = 18)
            ax2.set_ylabel('z ($10^3$ kilometers)', fontsize = 18)

            ax3.set_xlabel('z ($10^3$ kilometers)', fontsize = 18)
            ax3.set_ylabel('y ($10^3$ kilometers)', fontsize = 18)

            ax4.set_xlabel('x ($10^3$ kilometers)', fontsize = 18)
            ax4.set_ylabel('y ($10^3$ kilometers)', fontsize = 18)

            ax2.plot(nx,nz,label=names[1])
            ax4.plot(nx,ny,label=names[1])
            ax3.plot(nz,ny,label=names[1])

        else:
            ax2.plot(nx,nz,label=names[1])
            ax4.plot(nx,ny,label=names[1])
            ax3.plot(nz,ny,label=names[1])

        fig2.legend()
        fig2.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()
        pdf.savefig(fig2)
        
        plt.close(fig2)
        print("Keplerian figures saved to the local file as "+filename)
    return()
    
    
    
   