#
#	mm_synth.py
#
#	Generates synthetic astrometry to test multimoon
#
#	Benjamin Proudfoot
#	06/19/20
#

# I think the best method would be to have a tests directory where we store the synthetic astrometry,
# an output file associated with the creation of that synthetic astrometry (basically the parameters),
# and everything else needed to test the system (geocentric_obj_position.csv, etc).

# Starting psuedocode

# load in all required packages and functions
import numpy as np
.
.
.

# maybe have a runprops type thing to hold the inputs?

# change mm_chisquare to have an output option which just outputs the Model_DeltaLong and 
# Model_DeltaLat
# This will contain all of the data relevant to creating the synthetic astrometry

# Need to think about excluding data points where the secondary/tertiary/etc is aligned with the 
# primary. In theory the fitter will have no problem fitting to this data, but it really isn't 
# realistic...

# Once the model dfs are output, then all that needs to be done is to output it as a data file.

# I think the best way to do all of this is to use a REAL observations file to get the date/time
# of observations and the error bars, and then replace the data points with synthetic ones.
