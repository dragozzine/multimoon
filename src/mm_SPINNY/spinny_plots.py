import numpy as np
import time
from time import ctime
from quaternion import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

#TO DO: even if "don't save" is chosen, it saves a blank pdf
    
def spinny_plot(plot_df):    
    
    #N = int(((len(plot_df.columns.values)-3)/10))
    cols = plot_df.columns.values
    N = sum(1 for s in cols if 'X_Pos_' in s)

    t_arr = plot_df['Times'].values.flatten()
    
    names = np.empty(N,dtype='object')
    
    for n in range(0,N):
        name_x = plot_df.columns.values[1+n]
        names[n] = name_x.replace("X_Pos_","")

        globals()[names[n]+'_x'] = plot_df["X_Pos_"+names[n]]
        globals()[names[n]+'_y'] = plot_df["Y_Pos_"+names[n]]
        globals()[names[n]+'_z'] = plot_df["Z_Pos_"+names[n]]
        
        if globals()[names[n]+'_x'][0] == 0.0:
            name_prim = names[n]
        else:
            print("")
            
    while names[0] != "Sun":
        names = np.roll(names, 1)
        
    for n in range(1,N):
        globals()[names[n]+'_i'] = plot_df["inclination_"+names[n]]
        globals()[names[n]+'_e'] = plot_df["eccentricity_"+names[n]]
        globals()[names[n]+'_a'] = plot_df["semimajor_axis_"+names[n]]
        globals()[names[n]+'_O'] = plot_df["longitude_ascending_"+names[n]]
        globals()[names[n]+'_w'] = plot_df["argument_periapse_"+names[n]]
        globals()[names[n]+'_M'] = plot_df["mean_anomaly_"+names[n]]
        
        globals()[names[n]+'_obliq'] = plot_df["obliquity_"+names[n]]
        globals()[names[n]+'_prec'] = plot_df["precession_"+names[n]]
        globals()[names[n]+'_long'] = plot_df["longitude_"+names[n]]
        
        globals()[names[n]+'_spin'] = plot_df["spin_rate_"+names[n]]
        
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
    filename = name_prim+"_figures_"+t_current+".pdf"
    
    with PdfPages(filename) as pdf:
        for n in range(1,N):
            t = t_arr/(3600*24)
            if globals()[names[n]+'_x'][0] == 0.0:
                ax2.plot([0.0],[0.0], marker='o', color='black',label=names[n])
                ax3.plot([0.0],[0.0], marker='o', color='black',label=names[n])
                ax4.plot([0.0],[0.0], marker='o', color='black',label=names[n])
                fig2.suptitle(names[n]+' System Orbits',fontsize=18,fontweight='bold')
                
            elif globals()[names[n]+'_x'][0] != 0.0:
                fig1, ax1 = plt.subplots(3,2,sharex=True,figsize=(10,10))
                fig0, ax0 = plt.subplots(3,1,sharex=True,figsize=(10,10))
                
                if t[-1] > 1000:
                    t = t/365.25
                    ax1[2,0].set_xlabel('Time (years)')
                    ax1[2,1].set_xlabel('Time (years)')
                    ax0[2].set_xlabel('Time (years)')

                else:
                    ax1[2,0].set_xlabel('Time (days)')
                    ax1[2,1].set_xlabel('Time (days)')
                    ax0[2].set_xlabel('Time (days)')
                    
                ax1[0,0].plot(t,globals()[names[n]+'_i'],label=names[n])
                ax1[0,1].plot(t,globals()[names[n]+'_e'],label=names[n])
                ax1[1,0].plot(t,globals()[names[n]+'_O'],label=names[n])
                ax1[1,1].plot(t,globals()[names[n]+'_a'],label=names[n])
                ax1[2,0].plot(t,globals()[names[n]+'_M'],label=names[n])
                ax1[2,1].plot(t,globals()[names[n]+'_w'],label=names[n])
                
                ax0[0].plot(t,globals()[names[n]+'_obliq'],label=names[n])
                ax0[1].plot(t,globals()[names[n]+'_prec'],label=names[n])
                ax0[2].plot(t,globals()[names[n]+'_long'],label=names[n])

                nx = globals()[names[n]+'_x']
                ny = globals()[names[n]+'_y']
                nz = globals()[names[n]+'_z']
                
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
                    
                    ax2.plot(nx,nz,label=names[n])
                    ax4.plot(nx,ny,label=names[n])
                    ax3.plot(nz,ny,label=names[n])
                
                else:
                    ax2.plot(nx,nz,label=names[n])
                    ax4.plot(nx,ny,label=names[n])
                    ax3.plot(nz,ny,label=names[n])
     
                     
                ax1[0,0].set_title('Inclination')
                ax1[0,1].set_title('Eccentricity')
                ax1[1,0].set_title('Longitude of Ascending Node')
                ax1[1,1].set_title('Semi-major Axis')
                ax1[2,0].set_title('Mean Anomaly')
                ax1[2,1].set_title('Argument of Periapsis')

                ax1[0,0].set_ylabel('Degrees (from the ecliptic)')
                ax1[0,1].set_ylabel('')
                ax1[1,0].set_ylabel('Degrees')
                ax1[1,1].set_ylabel('Kilometers')
                ax1[2,0].set_ylabel('Degrees')
                ax1[2,1].set_ylabel('Degrees')
                
                ax1[0,0].grid()
                ax1[0,1].grid()
                ax1[1,0].grid()
                ax1[1,1].grid()
                ax1[2,0].grid()
                ax1[2,1].grid()
        
                fig1.suptitle('Orbital Parameters -- '+names[n],fontsize=18,fontweight='bold')
                fig1.tight_layout(rect=[0, 0.03, 1, 0.95])     
                
                ax0[0].set_title('Axial Obliquity')
                ax0[1].set_title('Axial Precession')
                ax0[2].set_title('Axial Longitude')#'Spin Rate')
                
                ax0[0].set_ylabel('Degrees')
                ax0[1].set_ylabel('Degrees')
                ax0[2].set_ylabel('Seconds$^{-1}$')
                
                ax0[0].ticklabel_format(useOffset=False,style='plain',axis='y')
                ax0[1].ticklabel_format(useOffset=False,style='plain',axis='y')
                ax0[2].ticklabel_format(useOffset=False,style='plain',axis='y')

                
                ax0[0].grid()
                ax0[1].grid()
                ax0[2].grid()
                
                fig0.suptitle('Spin Axis Orientation -- '+names[n],fontsize=18,fontweight='bold')
                fig0.tight_layout(rect=[0, 0.03, 1, 0.95]) 
  
                pdf.savefig(fig1)
                pdf.savefig(fig0)
                       
            # end of for loop
            
        fig2.legend()
        fig2.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()
        pdf.savefig(fig2)
        print("Figures saved to the local file as "+filename)
        
        plt.close(fig2)
       
    return()
             
        
        
        
        
        
        
        
        