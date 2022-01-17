# multimoon
Fit solar system small body multiple systems with non-Keplerian orbits including quadrupole shapes. 
Darin Ragozzine, Benjamin Proudfoot, Seth Pincock, Dallin Spencer
Brigham Young University
contact darin_ragozzine@byu.edu for more information

An extension of the code used in Ragozzine & Brown 2009 with many additions and improvements. Uses SPINNY (Simon Porter, Darin Ragozzine, et al.) as 
the n-quadrupole integrater at the core. Designed for Kuiper Belt Object systems (but ostensibly could be used for any solar system binary).

To install MultiMoon:

We highly recommend that you use a virtual environment to run MultiMoon. Once you have a virtual environment running, all package requirements are contained in `requirements.txt`
Next, inside of `multimoon/src/mm_SPINNY` execute `python setup.py build_ext --inplace`. This installs SPINNY in the appropriate location.

Talk about runprops.txt and all of the lines involved. (How to use)
