import commentjson as json
import pandas as pd
import sys
import os
import datetime
import shutil
import numpy as np

# DS TODO: Comment this

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
if 'runs' in cwd:
    runs_file = os.path.basename(os.path.normpath(cwd))

    os.chdir('../')
    cwd = os.getcwd()
    runs_file2 = os.path.basename(os.path.normpath(cwd))
    os.chdir('../../src')
    runprops['runs_file'] = '../runs/'+runs_file2+'/'+runs_file
    #print(runprops)

elif 'results' in cwd:
    runs_file = os.path.basename(os.path.normpath(cwd))
    os.chdir('../../../src')
    runprops['runs_file'] = '../results/'+objname+'/'+runs_file


if sys.argv[0] == "mm_synth_unseen.py" or sys.argv[0] == "mm_synth.py":
    runprops['first_run'] = False

if runprops.get("first_run") == True:
    import git
    repo = git.Repo(search_parent_directories=True)
    sha = repo.head.object.hexsha
    
    runprops['MultiMoon commit hash'] = sha
    
    x = datetime.datetime.now()

    date = str(x.strftime("%Y"))+"-"+str(x.strftime("%m"))+"-"+str(x.strftime("%d"))+"_"+str(x.strftime("%H"))+"."+str(x.strftime("%M"))+"."+str(x.strftime("%S"))
    runprops["date"] = date
    runprops['first_run'] = False
    
    
    now = datetime.datetime.now()

    time = now.strftime("%H:%M:%S")
    
    runprops["time"] = time
    

    newpath = "../results/"+objname+"/"+objname+"_"+date+"_"+runprops.get("run_file")
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    
    runprops['results_folder'] = newpath
    
    if isinstance(runprops.get('numpy_seed'), int):
        np.random.seed(seed=runprops.get('numpy_seed'))
    else:
        seed = x.microsecond*x.second		# This now creates a "random" seed based on the time
        np.random.seed(seed=seed)
        runprops['numpy_seed'] = seed
    
    if ('runs' in runprops.get('runs_file')): 
        shutil.copy(runprops.get('runs_file')+'/runprops.txt', newpath+'/runprops.txt')
    elif ('results' in runprops.get('runs_file')): 
        shutil.copy(runprops.get('runs_file')+'/runprops.txt', newpath+'/runprops.txt')
        shutil.copy(runprops.get('runs_file')+'/chain.h5', newpath+'/chain.h5')
        shutil.copy(runprops.get('runs_file')+'/best_likelihoods.csv', newpath+'/best_likelihoods.csv')
        runprops['chain_file'] = newpath+'/chain.h5'
    else:
        shutil.copy(filename, newpath+'/runprops.txt')
    
    init = runprops.get('runs_file')+'/'+objname+'_init_guess.csv'
    priors = runprops.get('runs_file')+'/'+objname+'_priors_df.csv'
    #geoanalysis = runprops.get('runs_file')+'/geocentric_'+objname+'_position_analysis.csv'

    if 'runs' in runprops.get('runs_file'):
        obs = runprops.get('runs_file')+'/../observations/'+runprops.get("obs_df")
    elif 'results' in runprops.get('runs_file'):
        obs = runprops.get('runs_file')+'/'+runprops.get('objectname')+'_obs_df.csv'

    runprops['init_filename'] = init
    runprops['obsdata_file'] = obs
    runprops['priors_filename'] = priors    
    
    shutil.copy(obs, newpath+'/'+runprops.get("objectname")+'_obs_df.csv')
    shutil.copy(priors, newpath+'/'+runprops.get("objectname")+'_priors_df.csv')
    shutil.copy(init, newpath+'/'+runprops.get("objectname")+'_init_guess.csv')
    #shutil.copy(geoanalysis, newpath+'/geocentric_'+runprops.get("objectname")+'_position_analysis.csv')
