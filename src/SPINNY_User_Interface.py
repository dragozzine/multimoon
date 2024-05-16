from mm_SPINNY.spinny_plots import spinny_plot
from mm_SPINNY.spinny_generate import *
from mm_SPINNY.spinny_nosun import *
from mm_SPINNY.mm_vpython import *
from mm_SPINNY.keplerian import *
from mm_SPINNY.spinny_vector import *
import numpy as np
import time
from time import ctime
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import sys
import commentjson as json

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

def main_menu():
    print("S: Run a system though SPINNY")
    print("F: Generate figures from an existing SPINNY integration")
    print("V: Generate a VPython animation")
    #print("C: Convert SPINNY data to/from barycentric and primaricentric frames") #### add later
    print("Q: Quit")
    
    user_input = str(input("Please make a selection: "))
    user_input = user_input.upper()
    
    if user_input == "S":
        run_spinny()
        
    elif user_input == "F":
        file = str(input("Please input SPINNY data filename (MUST be .csv): "))
        file = "output/"+file
        if not ".csv" in file:
            file += ".csv"
            
        plot_df = pd.read_csv(str(file),index_col=[0])
        plot(plot_df)
    elif user_input == "V":
        run_vpython()
            
    elif user_input == "Q":
        user_input = str(input("Are you sure you want to quit? (Y/N): "))
        user_input = user_input.upper()
        if user_input == "Y":
            print("Goodbye.")
            sys.exit()
        elif user_input == "N":
            print("\nReturning to main menu...")
            return main_menu()
        else:
            print("")
            print('Invalid Response.')
            return main_menu()

    else:
        print("")
        print('Invalid Response.')
        return main_menu()
        
                
def run_spinny():
    
    file_sys = str(input("Please input file name with system parameters (MUST be .csv): "))

    file_t = str(input("All times files have been moved to src/mm_SPINNY/times/. Enter name of file in that folder without the path. File name with list of times (MUST be .csv): "))
    
    #file_r = str(input("Input a runprops dictionary (MUST be .txt): "))
    print("Now using ../runs/runprops by default.")
    file_r = "../runs/runprops"
    
    if not ".csv" in file_sys:
        file_sys += ".csv"  
        
    if not ".csv" in file_t:
        file_t += ".csv"
        
    if not ".txt" in file_r:
        file_r += ".txt"
        
    #pd.set_option('display.max_columns', None)
    sys_df = pd.read_csv(str(file_sys),index_col=[0])
    time_df = pd.read_csv(str("mm_SPINNY/times/"+file_t))
    getData = ReadJson(file_r)
    runprops = getData.outProps()
    
    t_arr = time_df.values.flatten()
    N = len(sys_df.columns)
    
    j2_sum = sum(sys_df.loc["j2r2",:].values.flatten())
    names = list(sys_df.columns)
    
    tol = runprops.get("spinny_tolerance")
    
    if N == 2 and j2_sum == 0.00:
        kepler_system = kepler_integrate(sys_df,t_arr)
        kepler_df = kepler_system[0]
        names = kepler_system[1]
        kepler_save(kepler_df, names)
        print("\n Returning to main menu...")    
        return main_menu()
    elif not "Sun" in names:
        system = build_spinny_ns(sys_df,runprops)
        #system = build_spinny_multimoon(sys_df,runprops)
        spinny = evolve_spinny_ns(system[0],system[1],system[2],system[3],system[4],system[5],t_arr,tol,runprops)
        s_df = spinny[0]
        names = spinny[2]
        save(s_df,names, runprops)
    else:
        print('sys_df line 113', sys_df)
        system = build_spinny(sys_df,runprops)
        spinny = evolve_spinny(system[0],system[1],system[2],system[3],system[4],system[5],t_arr,runprops)
        s_df = spinny[0]
        print('s_df 117', s_df)
        names = spinny[1]
        save(s_df,names, runprops)

        
def save(s_df,names, runprops):

    save_yn = str(input("Do you want to save this data? (Y/N): "))
    save_yn = save_yn.upper()
    if save_yn=="Y":
        print("Generating .csv...")
        t_current = ctime().replace(" ","_")
        filename = names[1]+"_SPINNY_"+t_current+".csv"
        s_df.to_csv("output/"+filename)
        print("SPINNY data saved to the output file as "+filename)
        plot_q(s_df, names, runprops)
    elif save_yn == "N":
        print("")
        plot_q(s_df, names, runprops)
    else:
        print("")
        print('Invalid Response.')
        return save(s_df,names)   
    
    
def plot_q(s_df, names, runprops): 
    plot_yn = str(input("Do you want me to generate figures from these data? (Y/N): "))
    plot_yn = plot_yn.upper()
    
    if plot_yn == "Y":
        plot(s_df, names, runprops)
    elif plot_yn == "N":
        print("\n Returning to main menu...")
        return main_menu()
    else:
        print("")
        print('Invalid Response.')
        return plot_q(s_df, names, runprops)
           

def plot(plot_df, names, runprops):

    print("\n Generating figures...")
    spinny_plot(plot_df, names, runprops)
    print("\n Returning to main menu...")
    return main_menu()


def run_vpython():
    file = str(input("Please input file name with parameters of system to animate (MUST be .csv): "))
    if not ".csv" in file:
        file += ".csv"  

    sys_df = pd.read_csv(str(file),index_col=[0])
    print("Which data do you want me to animate?")
    print("S: Spins")
    print("O: Orbits")
    user_input = str(input("Please select : "))
    
    user_input = user_input.upper()
    if user_input == "O":
        vpython_orbits(sys_df)
        return main_menu()
    elif user_input == "S":
        vpython_spins(sys_df)
        print("\nReturning to main menu...")    
        return main_menu()
    else:
        print("")
        print('Invalid Response.')
        return run_vpython()
    
print("Hello!")
print("I am the SPINNY user interface. What can I help you with?")
main_menu()





















