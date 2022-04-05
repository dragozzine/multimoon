
import pandas as pd
import os
from tqdm import tqdm
from latlon_transform import *
from mm_make_geo_pos import *
import sys
import shutil
import commentjson as json
from scipy.spatial.transform import Rotation as R

class ReadJson(object):
    def __init__(self, filename):
        self.data = json.load(open(filename))
    def outProps(self):
        return self.data


