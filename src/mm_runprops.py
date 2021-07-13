import commentjson as json
import pandas as pd
import sys
import os
import datetime
import shutil

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data
runs_file = ''
cwd = ''
filename = ""
cwd = os.getcwd()
if len(sys.argv) > 1:
    filename = sys.argv[1]
elif sys.argv[0] == "mm_synth.py":
    filename = "runprops_gensynth.txt"
elif sys.argv[0] == "mm_synth_unseen.py":
    filename = "runprops_gensynth_unseen.txt"
else:
    #cwd = os.getcwd()
    if 'src' in cwd:
        filename = "../runs/runprops.txt"
    elif 'runs' in cwd:
        filename = "runprops.txt"
    elif 'results' in cwd:
        filename = "runprops.txt"
        
    else:
        print('You are not starting from a proper directory, You should run mm_run.py from either a runs directory or a results directory.')
        sys.exit()
    
getData = ReadJson(filename)
runprops = getData.outProps()
runprops['chain_file'] = None
runprops['first_run'] = True
objname = runprops.get("objectname")
#print(cwd)
if 'runs' in cwd:
    runs_file = os.path.basename(os.path.normpath(cwd))
    os.chdir('../../../src')
    runprops['runs_file'] = '../runs/'+objname+'/'+runs_file
    #print(runprops)
elif 'results' in cwd:
    runs_file = os.path.basename(os.path.normpath(cwd))
    os.chdir('../../../src')
    runprops['runs_file'] = '../results/'+objname+'/'+runs_file


if sys.argv[0] == "mm_synth_unseen.py" or sys.argv[0] == "mm_synth.py":
    runprops['first_run'] = False

if runprops.get("first_run") == True:
    x = datetime.datetime.now()

    date = str(x.strftime("%Y"))+"-"+str(x.strftime("%m"))+"-"+str(x.strftime("%d"))+"_"+str(x.strftime("%H"))+"."+str(x.strftime("%M"))+"."+str(x.strftime("%S"))
    runprops["date"] = date
    runprops['first_run'] = False
    

    newpath = "../results/"+objname+"/"+objname+"_"+date+"_"+runprops.get("run_file")
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    
    runprops['results_folder'] = newpath
    
    #print('runs file: ', runprops.get('runs_file'))
    #print('results file: ', runprops.get('results_folder'))
    
    if ('runs' in runprops.get('runs_file')): 
        shutil.copy(runprops.get('runs_file')+'/runprops.txt', newpath+'/runprops.txt')
    elif ('results' in runprops.get('runs_file')): 
        shutil.copy(runprops.get('runs_file')+'/runprops.txt', newpath+'/runprops.txt')
        shutil.copy(runprops.get('runs_file')+'/chain.h5', newpath+'/chain.h5')
        runprops['chain_file'] = newpath+'/chain.h5'
    else:
        shutil.copy(filename, newpath+'/runprops.txt')
    #shutil.copy(filename, '../runs/'+objname+'/runprops.txt')
    
    init = runprops.get('runs_file')+'/'+objname+'_init_guess.csv'
    priors = runprops.get('runs_file')+'/'+objname+'_priors_df.csv'
    geoanalysis = runprops.get('runs_file')+'/geocentric_'+objname+'_position_analysis.csv'
    print("\n\nWarning: Observations data frames are now centrally located. Loading the centrally located obs_df from runs/"+objname+"/observations/.")
    #print(sys.cwd())

    if 'runs' in runprops.get('runs_file'):
        obs = runprops.get('runs_file')+'/../observations/'+runprops.get("obs_df")
    #obs = runprops.get('runs_file')+'../observations/'+objname+'_obs_df.csv'
    elif 'results' in runprops.get('runs_file'):
        obs = runprops.get('runs_file')+'/'+runprops.get('objectname')+'_obs_df.csv'

    #print(init,priors,obs)
    
    runprops['init_filename'] = init
    runprops['obsdata_file'] = obs
    runprops['priors_filename'] = priors    
    
    shutil.copy(obs, newpath+'/'+runprops.get("objectname")+'_obs_df.csv')
    shutil.copy(priors, newpath+'/'+runprops.get("objectname")+'_priors_df.csv')
    shutil.copy(init, newpath+'/'+runprops.get("objectname")+'_init_guess.csv')
    shutil.copy(geoanalysis, newpath+'/geocentric_'+runprops.get("objectname")+'_position_analysis.csv')
    #runprops["init_filename"] = "../data/" + runprops.get("objectname") + "/" + runprops.get("objectname") + "_init_guess.csv"
    #runprops["priors_filename"] = "../data/" + runprops.get("objectname") + "/" + runprops.get("objectname") + "_priors_df.csv"
    #runprops["obsdata_file"] = "../data/" + runprops.get("objectname") + "/" + runprops.get("objectname") + "_obs_df.csv"
