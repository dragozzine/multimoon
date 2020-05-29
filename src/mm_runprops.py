import json
import pandas as pd

class ReadJson(object):
    def __init__(self):
        print('Read the runprops.txt file')
        self.data = json.load(open("runprops.txt"))
    def outProps(self):
        return self.data

getData = ReadJson()
runprops = getData.outProps()

runprops["init_filename"] = "../data/" + runprops.get("objectname") + "/" + runprops.get("objectname") + "_init_guess.csv"
runprops["priors_filename"] = "../data/" + runprops.get("objectname") + "/" + runprops.get("objectname") + "_priors_df.csv"
runprops["obsdata_file"] = "../data/" + runprops.get("objectname") + "/" + runprops.get("objectname") + "_obs_df.csv"

fixfloat_df = pd.DataFrame(data = runprops.get('float_dict'), index = [0])

lpriorfunc = 'log_prior'
llhoodfunc = 'log_likelihood'

posteriorfile = 'date_time_objname.csv'
verbose = 1
plotflags = True
plot_model_in_data = False
plot_model_over_time = False
