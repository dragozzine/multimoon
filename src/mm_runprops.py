import json
import pandas as pd
import sys


class ReadJson(object):
    def __init__(self, filename):
        print('Read the runprops.txt file')
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data

if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = "runprops.txt"
    
getData = ReadJson("runprops.txt")
runprops = getData.outProps()

runprops["init_filename"] = "../data/" + runprops.get("objectname") + "/" + runprops.get("objectname") + "/_init_guess.csv"
runprops["priors_filename"] = "../data/" + runprops.get("objectname") + "/" + runprops.get("objectname") + "/_priors.csv"
runprops["obsdata_file"] = "../data/" + runprops.get("objectname") + "/" + runprops.get("objectname") + "/ObsDF.csv"