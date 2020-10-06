from astroquery.jplhorizons import Horizons
import numpy as np
import pandas as pd
from astropy.time import Time
import time

def mm_make_geo_pos(objname, times, runprops, synthetic = False):
    """This function takes a name of a solar system body(a KBO), and creates a csv file of the body's ephemerides"""
    starttime = time.time()
    
    ourKBO = Horizons(id=objname,location=399,epochs = times)
    ephKBO = ourKBO.ephemerides()['RA','DEC','datetime_jd']
    vecKBO = ourKBO.vectors()['lighttime','x','y','z']
    
    jdTime= ephKBO['datetime_jd']
    lightTime = vecKBO['lighttime']
    kboTime=jdTime-lightTime
    
    #taking care to convert the AU distances to kilometers
    geocentricFile = pd.DataFrame({'geocentricTime':ephKBO['datetime_jd'],'x':vecKBO['x']*149597870.7,'y':vecKBO['y']*149597870.7 ,'z':vecKBO['z']*149597870.7})
    outFile = pd.DataFrame({'kboTIME':kboTime,'x':vecKBO['x']*149597870.7 ,'y':vecKBO['y']*149597870.7 ,'z':vecKBO['z']*149597870.7 })

    if synthetic:
        return outFile

    filename1 = "../runs/"+objname+"/"+runprops.get('run_file')+"/geocentric_" + objname + "_position.csv"
    geocentricFile.to_csv(filename1)
    #fileName2 = "geocentric_" + objectName +'_position.csv'
    #outFile.to_csv(fileName2)
    
    endtime = time.time()
    
    runLength = endtime-starttime
    #print("Runtime was: "+str(runLength)+" seconds" )
    
'''
    NAME:
         geotoKBOtime
         
    PURPOSE:
         Given the name of an object, this function will open a file bearing that object's name,
         and converts the julian data in that file to Heliocentric JD.
         
    CALLING SEQUENCE:
         geotoKBOtime(objectName)
   
    INPUTS
          objectName - The name of the KBO (or Horizons id). The name has 2 rules.
          1) The name is searchable in JPL Horizons
          2) The name is the beginning of a file titled "Data_geoTime.csv" in the same directory
          
    OUTPUTS:
          No functions outputs, but does create a new csv file for the object given.
'''
def geotoKBOtime(objectName):
    starttime = time.time()
    theFile = "../data/" + objectName + "/" + objectName+'Data_geoTime.csv'
    
    oldData = pd.read_csv(theFile)
    date = oldData['Dates']
    
    date1 = date.iloc[0]
    date2 = date.iloc[-1]
    
    jddate1 = Time(date1,format='jd')
    jddate2 = Time(date2,format='jd')
    
    dateList = []
    
    for i in date:
        dateList.append(Time(i,format='jd'))
                        
    theObject = Horizons(id=objectName, location = 399,epochs=dateList)
    vecObj = theObject.vectors()['lighttime']
    
    kboTime=date-vecObj
    
    
    newData = oldData
    newData['Dates'] = kboTime
    
    newData=newData.rename(columns = {'Dates':'datetime_jd'})
    
    outFile = "../data/" + objectName + "/"+ objectName+'Data_KBOTime.csv'
   
    newData.to_csv(outFile)
    
    endtime = time.time()
    
    runLength = endtime-starttime
    print("Runtime was: "+str(runLength)+" seconds" )
    
