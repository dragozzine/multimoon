import commentjson as json
import pandas as pd
import sys
import os
import datetime

class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

filename = ""
if len(sys.argv) > 1:
    filename = sys.argv[1]
elif sys.argv[0] == "mm_synth.py":
    filename = "runprops_gensynth.txt"
else:
    filename = "../runs/runprops.txt"
    
getData = ReadJson(filename)
runprops = getData.outProps()

if runprops.get("first_run") == True:
    x = datetime.datetime.now()

    date = str(x.strftime("%Y"))+"-"+str(x.strftime("%m"))+"-"+str(x.strftime("%d"))+"_"+str(x.strftime("%H"))+"."+str(x.strftime("%M"))+"."+str(x.strftime("%S"))
    runprops["date"] = date
    runprops['first_run'] = False
    
    newpath = "../runs/"+runprops.get("objectname")+"_"+date
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    
    runprops['runs_folder'] = newpath
    import shutil
    shutil.copy(filename, '../runs/'+newpath+'/runprops.txt')
    
    init = runprops.get('init_filename')
    priors = runprops.get('priors_filename')
    obs = runprops.get('obsdata_file')
    
    shutil.copy(obs, '../runs/'+newpath+'/'+runprops.get("objectname")+'_obs_df.csv')
    shutil.copy(priors, '../runs/'+newpath+'/'+runprops.get("objectname")+'_priors_df.csv')
    shutil.copy(init, '../runs/'+newpath+'/'+runprops.get("objectname")+'_init_guess.csv')
    
    #runprops["init_filename"] = "../data/" + runprops.get("objectname") + "/" + runprops.get("objectname") + "_init_guess.csv"
    #runprops["priors_filename"] = "../data/" + runprops.get("objectname") + "/" + runprops.get("objectname") + "_priors_df.csv"
    #runprops["obsdata_file"] = "../data/" + runprops.get("objectname") + "/" + runprops.get("objectname") + "_obs_df.csv"
