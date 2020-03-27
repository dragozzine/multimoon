# multimoon
Fit solar system small body multiple systems with non-Keplerian orbits including quadrupole shapes. 
Darin Ragozzine, Benjamin Proudfoot, Seth Pincock, Dallin Spencer
Brightam Young University
contact darin_ragozzine@byu.edu for more information

An extension of the code used in Ragozzine & Brown 2009 with many additions and improvements. Uses SPINNY (Simon Porter, Darin Ragozzine, et al.) as the n-quadrupole integrater at the core. Designed for Kuiper Belt Object systems (but ostensibly could be used for any solar system binary). 

FOR NOW: 
to run these codes use the following at the top of the code:

import sys
sys.append.path('./src')

or whatever directory your MultiMoon "src" directory is

then you can import individual files and functions. Example: 

import testsupport
out=testsupport.test_runprops()
