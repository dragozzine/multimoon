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

print(runprops.get("parameter_dict"))

runprops.get('priors_filename') = "../data/" + runprops.get("objectname") + "/" + runprops.get("objectname") + "/_priors.csv"
runprops.get("obsdata_file") = "../data/" + runprops.get("objectname") + "/" + runprops.get("objectname") + "/ObsDF.csv"


float_dict = {'name_1': 0,
              'mass_1': 0,
              'name_2': 0,
              'mass_2': 0,
              'sma_2': 0,
              'ecc_2': 0,
              'aop_2': 0,
              'inc_2': 0,
              'lan_2': 0,
              'mea_2': 0,
              'name_3': 0,
              'mass_3': 0,
              'sma_3': 0,
              'ecc_3': 0,
              'aop_3': 0,
              'inc_3': 0,
              'lan_3': 0,
              'mea_3': 0,
              'j2r2_1': 0,
              'c22r2_1': 0,
              'spaop_1': 0,
              'spinc_1': 0,
              'splan_1': 0,
              'sprate_1': 0}
fixfloat_df = pd.DataFrame(data = float_dict, index = [0])

lpriorfunc = 'log_prior'
llhoodfunc = 'log_likelihood'

posteriorfile = 'date_time_objname.csv'
verbose = 1
plotflags = True
plot_model_in_data = False
plot_model_over_time = False