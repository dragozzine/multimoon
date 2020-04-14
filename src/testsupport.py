# testsupport.py
#
# Darin Ragozzine
# March 27, 2020
# 
# Supports testing of MultiMoon by including functions
# that return valid "shapes" and types of various
# planned pieces of MultiMoon. 

import pandas
from astroquery.jplhorizons import Horizons

# test_runprops
# 

def test_runprops():
    """ Returns a runprops dictionary. """
    
    runprops = {
        'date': 'today',
        'time': 'now',
        'user': 'autogenerate',
        'MultiMoon commit hash': 'autogenerate',
        'mm_initguess_which' : 'get_init_guess_from_dist', # which function from mm_initguess to use
        'numobjects' : 3,
        'RunGoal' : 'user puts their goal for this run here',
        'objectname' : 'Haumea',
        'nwalkers' : 30,
        'nsteps' : 1001,
        'nthinning' : 100
    }
    
#Filename for starting guess - default = start_guess.csv
#fixfloat_df - default="all float" dataframe in parameters format
#lpriorfunc - default="log_prior"
#llhoodfunc - default="log_likelihood"
#obsdata - name of file with text/csv of observations dataframe
#posteriorfile - name of file for posterior, default="date_time_objname.csv"
#verbose - default=1, how verbose to be in terms of output
#(0 = quiet, no text; 1 = typical, only major outputs, 5 = everything) 
#plotflags: True if this plot is to be made, False if not
##plot_model_in_data - show model in data space
#plot_model_over_time - show model continuous in time over data space
#other plots too

    return runprops


def test_mmparamdf():
    """ Returns a test parameters dataframe with fake values."""

    objects=["0","1","2","3"] # 3 objects in KBO system and Sun (0)
    objparams=["name","mass","sma","ecc","aop","inc","lan","mea","j2r2","c22r2","spaop","spinc","splan","sprate"]
    
    columns=[]
    
    for objectnum in objects:
        for param in objparams:
            columns.append(param+"_"+objectnum)

    mmparamdf=pandas.DataFrame(index=[0],columns=columns)
    
    mmparamdf["name_0"]="Sun"
    mmparamdf["mass_0"]=1.0
    
    
    mmparamdf.fillna(0.0)
    
    return mmparamdf


def mm_make_geo_pos(objname, start='2000-01-01', end='2040-01-01', step='10d'):
    obj = Horizons(id=objname,location='500',epochs={'start':start, 'stop':end,'step':step})
    obsDF = pandas.DataFrame()
    obsDF['time'] = obj.vectors()['datetime_jd']
    obsDF['x'] = obj.vectors()['x']
    obsDF['y'] = obj.vectors()['y']
    obsDF['z'] = obj.vectors()['z']
    
    filename = 'geocentric_'+objname+'_position.csv'
    obsDF.to_csv(filename,index=False)  


#test_runprops()
print(test_mmparamdf())
