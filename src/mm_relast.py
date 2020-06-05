import numpy as np
from astropy import units as u
import pandas as pd
from math import radians, sin, cos, sqrt, asin

'''
    NAME:
         convert_ecl_rel_pos_to_geo_rel_ast
         
    PURPOSE:
         To determine the deltaLat and deltaLong of a moon from its primary KBO.
         
    CALLING SEQUENCE:
          deltaLong,deltaLat = convert_ecl_rel_pos_to_geo_rel_ast(obs_J2000_pos, obj_J2000_pos, moon_obs_pos)
   
    INPUTS
          obs_J2000_pos - The J2000 ecliptic relative position of the Observer in Cartesian coordinates
          obj_J2000_pos - The J2000 ecliptic relative position of the KBO in Cartesian coordinates
          moon_obs_pos -  The observer relative position of the Moon in Cartesian coordinates
   
    OUTPUTS:
          deltaLong - The difference in Longitude of the moon vs. it's primary KBO in arcseconds
          deltaLat - The difference in Latitude of the moon vs. it's primary KBO in arcseconds
'''
def convert_J2000_ecl_to_geo_pos_rel(obs_J2000_pos, obj_J2000_pos, moon_obs_pos):
    
    #Get the Cartesian positions of the Observer
    x1,y1,z1 = ecl_rel_pos[0],ecl_rel_pos[1],ecl_rel_pos[2]
    
    #Get the distance from the origin (Heliocenter) of the observer
    R1= np.sqrt(x1**2+y1**2+z1**2)    

    #Get the Heliocentric Cartesian positions of the KBO
    x2,y2,z2 = obj_rel_pos[0],obj_rel_pos[1],obj_rel_pos[2]
    
    #Observer centric coordinates on KBO
    x2,y2,z2 = x2-x1,y2-y1,z2-z1
    
    #Get the distance from the origin (Now observer-center) of the KBO
    R2= np.sqrt(x2**2+y2**2+z2**2)
    
    x3,y3,z3 = rel_moon[0],rel_moon[1],rel_moon[2]
    moonX = x3+x2
    moonY = y3+y2
    moonZ = z3+z2
    R3 = np.sqrt((x3+x2)**2+(y3+y2)**2+(z3+z2)**2)
    
    
    #Now calculate the latitude and longitude from the coordinates given
    longitude1 = np.arcsin(z1/R1)*360/2/np.pi
    latitude1 = np.arccos(x1/R1/np.cos(longitude1*u.degree))/u.rad*360/2/np.pi
    
    #Calculate the latitude and longitude from the coordinates
    longitude2 = np.arcsin(z2/R2)*360/2/np.pi
    latitude2 = np.arccos(x2/R2/np.cos(longitude2*u.degree))/u.rad*360/2/np.pi
    
     #Calculate the latitude and longitude from the coordinates
    longitude3 = np.arcsin(moonZ/R3)*360/2/np.pi
    latitude3 = np.arccos(moonX/R3/np.cos(longitude3*u.degree))/u.rad*360/2/np.pi
    
    #Calculate the deltaLat and deltaLong
    deltaLat = (latitude3-latitude2)
    deltaLong = (longitude3-longitude2)*np.cos(latitude2*u.degree)

    return deltaLong*3600, deltaLat*3600

def convert_ecl_rel_pos_to_geo_rel_ast(obs_to_prim_pos, prim_to_sat_pos):
    
     #Get the Cartesian positions of the Observer
    x1,y1,z1 = obs_to_prim_pos[0],obs_to_prim_pos[1],obs_to_prim_pos[2]
    
    #Get the distance from the observer to the primary
    R1= np.sqrt(x1**2+y1**2+z1**2)    
    
    x2,y2,z2 = prim_to_sat_pos[0],prim_to_sat_pos[1],prim_to_sat_pos[2]
    
    x2,y2,z2 = x1+x2,y1+y2,z1+z2
    
    R2 = np.sqrt(x2**2+y2**2+z2**2)
    

    longitude1 = np.arcsin(z1/R1)*360/2/np.pi
    latitude1 = np.arccos(x1/R1/np.cos(longitude1*u.degree))/u.rad*360/2/np.pi

    longitude2 = np.arcsin(z2/R2)*360/2/np.pi
    latitude2 = np.arccos(x2/R2/np.cos(longitude2*u.degree))/u.rad*360/2/np.pi
    
    deltaLat = (latitude2-latitude1)
    deltaLong = (longitude2-longitude1)*np.cos(latitude1*u.degree)
    
    #ang_dist = haversine(latitude1, longitude1, latitude2, longitude2)
    
    return deltaLong*3600, deltaLat*3600

def haversine(lat1, lon1, lat2, lon2):
 
    dLat = radians(lat2 - lat1)
    dLon = radians(lon2 - lon1)
    lat1 = radians(lat1)
    lat2 = radians(lat2)
 
    a = sin(dLat / 2)**2 + cos(lat1) * cos(lat2) * sin(dLon / 2)**2
    c = 2 * asin(sqrt(a))
 
    return c
    