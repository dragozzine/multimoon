import numpy as np
import time
from time import ctime
from quaternion import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

    
def spinny_plot(plot_df, names):    
    
    cols = plot_df.columns.values
    N = sum(1 for s in cols if 'X_Pos_' in s)

    t_arr = plot_df['Times'].values.flatten()
    T = len(t_arr)
    names = names
    
    for name in names:

        globals()[name+'_x'] = plot_df["X_Pos_"+name]
        globals()[name+'_y'] = plot_df["Y_Pos_"+name]
        globals()[name+'_z'] = plot_df["Z_Pos_"+name]
        globals()[name+'_xL'] = plot_df["Lx_"+name]
        globals()[name+'_yL'] = plot_df["Ly_"+name]
        globals()[name+'_zL'] = plot_df["Lz_"+name]
        globals()[name+'_E'] = plot_df["E_"+name]
        
        if globals()[name+'_x'][0] == 0.0:
            name_prim = name
        else:
            print("")
            
        globals()[name+'_i'] = plot_df["inclination_"+name]
        globals()[name+'_e'] = plot_df["eccentricity_"+name]
        globals()[name+'_a'] = plot_df["semimajor_axis_"+name]
        globals()[name+'_O'] = plot_df["longitude_ascending_"+name]
        globals()[name+'_w'] = plot_df["argument_periapse_"+name]
        globals()[name+'_M'] = plot_df["mean_anomaly_"+name]
        
        globals()[name+'_obliq'] = plot_df["obliquity_"+name]
        globals()[name+'_prec'] = plot_df["precession_"+name]
        globals()[name+'_long'] = plot_df["longitude_"+name]
        
        globals()[name+'_spin_angle'] = plot_df["spin_orbit_angle_"+name]
        globals()[name+'_spin_rate'] = plot_df["spin_rate_"+name]
        
    
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
    filename = "../results/SPINNY-models/"+name_prim+"_figures_"+t_current+".pdf"
    
    with PdfPages(filename) as pdf:
        
        for name in names:
            t = t_arr/(3600*24)
            if name == name_prim: # plots the primary as only a black dot at the origin
                ax2.plot([0.0],[0.0], marker='o', color='black', label = name)
                ax3.plot([0.0],[0.0], marker='o', color='black')
                ax4.plot([0.0],[0.0], marker='o', color='black')
                fig2.suptitle(name+' System Orbits',fontsize=18,fontweight='bold')
                
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
            
        ##### PLOTS THE ORBITAL PARAMETERS #####
            if name != name_prim: # exclude the primary
                ax1[0,0].plot(t,globals()[name+'_i'],label=name)
                ax1[0,1].plot(t,globals()[name+'_e'],label=name)
                ax1[1,0].plot(t,globals()[name+'_O'],label=name)
                ax1[1,1].plot(t,globals()[name+'_a'],label=name)
                ax1[2,0].plot(t,globals()[name+'_M'],label=name)
                ax1[2,1].plot(t,globals()[name+'_w'],label=name)

                # this generates line of best fit for energy/momentum for calculation of fractional change, if needed
                trend_x = np.polyfit(t,globals()[name+'_xL'], 1)
                px = np.poly1d(trend_x)
                trend_y = np.polyfit(t,globals()[name+'_yL'], 1)
                py = np.poly1d(trend_y)
                trend_z = np.polyfit(t,globals()[name+'_zL'], 1)
                pz = np.poly1d(trend_z)

                trend_E = np.polyfit(t,globals()[name+'_E'], 1)
                pE = np.poly1d(trend_E)

                #plots trendline on energy/momentum plots if needed  
                #ax0[1].plot(t,px(t),"r--",t,py(t),"r--",t,pz(t),"r--")
                #ax0[2].plot(t,pE(t),"r--")

                L_tot = [np.linalg.norm([px(t)[i],py(t)[i],pz(t)[i]]) for i in range(0,T)]

                nx = globals()[name+'_x']
                ny = globals()[name+'_y']
                nz = globals()[name+'_z']

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

                    ax2.plot(nx,nz,label=name)
                    ax4.plot(nx,ny,label=name)
                    ax3.plot(nz,ny,label=name)

                else:
                    ax2.plot(nx,nz,label=name)
                    ax4.plot(nx,ny)
                    ax3.plot(nz,ny)


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

                ax1[0,0].ticklabel_format(useOffset=False,style='plain',axis='y')
                ax1[0,1].ticklabel_format(useOffset=False,style='plain',axis='y')
                ax1[1,0].ticklabel_format(useOffset=False,style='plain',axis='y')
                ax1[1,1].ticklabel_format(useOffset=False,style='plain',axis='y')
                ax1[2,0].ticklabel_format(useOffset=False,style='plain',axis='y')
                ax1[2,1].ticklabel_format(useOffset=False,style='plain',axis='y')

                fig1.suptitle('Orbital Parameters -- '+name,fontsize=18,fontweight='bold')
                fig1.tight_layout(rect=[0, 0.03, 1, 0.95]) 
                pdf.savefig(fig1)
            
        ##### PLOTS SPIN PARAMETERS #####
        # all bodies included #
        
            ax0[0].plot(t,globals()[name+'_obliq'],label=name)
            ax0[1].plot(t,globals()[name+'_prec'],label=name)
            ax0[2].plot(t,globals()[name+'_spin_rate'],label=name)
            
            ax0[0].set_title('Axial Obliquity')
            ax0[1].set_title('Axial Precession')
            ax0[2].set_title('Spin Rate')

            ax0[0].set_ylabel('Degrees')
            ax0[1].set_ylabel('Degrees')
            ax0[2].set_ylabel('Seconds$^{-1}$')

            ax0[0].ticklabel_format(useOffset=False,style='plain',axis='y')
            ax0[1].ticklabel_format(useOffset=False,style='plain',axis='y')
            ax0[2].ticklabel_format(useOffset=False,style='plain',axis='y')

            ax0[0].grid()
            ax0[1].grid()
            ax0[2].grid()

            fig0.suptitle('Spin Parameters -- '+name,fontsize=18,fontweight='bold')
            fig0.tight_layout(rect=[0, 0.03, 1, 0.95]) 
            
            pdf.savefig(fig0)
                       
            # end of for loop
            
        percent_changeL = (L_tot[-1]-L_tot[0])/L_tot[0]
        percent_changeE = (pE(t)[-1]-pE(t)[0])/pE(t)[0]
        print("Fractional change in E: "+str(percent_changeE))
        print("Fractional change in L: "+str(percent_changeL))
        
        ax2.legend()    
        fig2.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()
        pdf.savefig(fig2)
        print("Figures saved to resfile as "+filename)
        
        plt.close(fig2)
       
    return()
             
        
        
        
        
        
        
        
        