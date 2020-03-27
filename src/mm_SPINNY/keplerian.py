import numpy as np
import spiceypy as spice
import pandas as pd
from spinny_generate import vec2orb
from time import ctime
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D

# computes a two-body, Keplerian integration using spiceypy
def kepler_integrate(sys_df,t_arr):
    
    G = 6.67e-20 # Gravitational constant in km
    N = 2
    T = len(t_arr)
    
    names = [sys_df.columns[0],sys_df.columns[1]]
    masses = sys_df.loc["mass",names[1]] + sys_df.loc["mass",names[0]]   

    mu = G*masses # compute the gravitational parameter for the system
    
    a0 = sys_df.loc["sma",names[1]]
    e0 = sys_df.loc["ecc",names[1]]
    i0 = sys_df.loc["inc",names[1]]
    O0 = sys_df.loc["lan",names[1]]
    w0 = sys_df.loc["aop",names[1]]
    M0 = sys_df.loc["mea",names[1]]
    p0 = a0*(1-e0) # periapsis distance
    
    T0 = 0.0 # epoch (seconds past J2000)
    
    orb0 = [p0,e0,i0,O0,w0,M0,T0,mu]
    
    init_vec = spice.conics(orb0,0.0)
    
    vec_arr = np.empty((T,6))
    orb_arr = np.empty((T,8))
    
    start_time = time.time()
    print("Running two-body Keplerian integration...")
    for t in range(0,T):
    
        vec_arr[t] = spice.conics(orb0,t_arr[t])
        orb_arr[t] = spice.oscelt(vec_arr[t],t,mu)
                         
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
    
    body_dict.setdefault('X_Pos_'+names[1] , vec_arr[:,0])
    body_dict.setdefault('Y_Pos_'+names[1] , vec_arr[:,1])
    body_dict.setdefault('Z_Pos_'+names[1] , vec_arr[:,2])
    body_dict.setdefault('X_Vel_'+names[1] , vec_arr[:,3])
    body_dict.setdefault('Y_Vel_'+names[1] , vec_arr[:,4])
    body_dict.setdefault('Z_Vel_'+names[1] , vec_arr[:,5])
    
    body_dict.setdefault('sma_'+names[1] , orb_arr[:,0]/(1-orb_arr[:,1]))
    body_dict.setdefault('ecc_'+names[1] , orb_arr[:,1])
    body_dict.setdefault('inc_'+names[1] , orb_arr[:,2])
    body_dict.setdefault('lan_'+names[1] , orb_arr[:,3])
    body_dict.setdefault('aop_'+names[1] , orb_arr[:,4])
    body_dict.setdefault('mea_'+names[1] , orb_arr[:,5])

    kepler_df = pd.DataFrame(body_dict)
    
    save(kepler_df, names)
    

def save(kepler_df, names):   
    save_yn = str(raw_input("Do you want to save this data? (Y/N): "))
    save_yn = save_yn.upper()
    
    if save_yn=="Y":
        print("Generating .csv...")
        t_current = ctime().replace(" ","_")
        filename = names[0]+"_KEPLER_"+t_current+".csv"
        kepler_df.to_csv(filename)
        print("Data from SpiceyPy saved to the local file as "+filename)
        plot_q(kepler_df,names)

    elif save_yn == "N":
        print("")
        plot_q(kepler_df,names)
    else:
        print("")
        print('Invalid Response.')
        return save(kepler_df,names) 
    
def plot_q(kepler_df,names): 
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
    
    
    
   