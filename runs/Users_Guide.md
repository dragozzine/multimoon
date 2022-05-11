# A User's Guide to Using MultiMoon

### Purpose of MultiMoon

MultiMoon is a program designed to fit solar system small body multiple systems with non-Keplerian orbits including their quadrupole shapes. This will allow us to learn more about distant solar system objects, and determine their orbital parameters through Bayesian statistical sampling. 
This code requires several types of files to be able to run. In particular, to get MultiMoon to function properly, you will need to provide a...
1) A csv file of the relative astrometry of the object
2) A csv file with the position of the object in geocentric cartesian coordinates
3) A csv file of the inital guess for the physical and orbital parameters of the object
4) A csv file that contains the priors for the object
5) A csv file containing information about the run


In the `runs/Lempo` folder, we have two examples of complete runs directories which can work out of the box. 

### Observation CSV,  "Object"_obs_df.csv

The relative astrometry of the binary in question is contained in this file. As MultiMoon uses J2000 geocentric ecliptic coordinates in its fitting, astrometry should be provided in J2000 geocentric ecliptic coordinates. A script for conversion of relative astrometry in J2000 geocentric equatorial coordinates to the appropriate coordinates is located at `src/latlon_transform.py`. 

Generally, the obs_df files use the following format (with Lempo as an example):

###### time, DeltaLat_Hiisi, DeltaLong_Hiisi, DeltaLat_Hiisi_err, DeltaLong_Hiisi_err, DeltaLat_Paha, DeltaLong_Paha, DeltaLat_Paha_err, DeltaLong_Paha_err

In reality, the order of the columns does not matter, since the columns are called through the use of pandas dataframes. The names in the columns must match the names given in the `runprops.txt` file. 

Observations for objects at certain times can be left empty if data does not exist at that time (only really applicable for 3 body or higher systems). The observation dataframe that you would like to use must be specified in the runprops.txt file as the obsdata_file variable.

### Initial Guess, "Object"_init_guess.csv

The initial guess file allows the user to specify a distribution from which initial walker position will be drawn. For each parameter, two numbers must be supplied, the mean and standard deviation of the walker distribution. The user is responsible for ensuring that the correct parameters are used in the initial guess file; MultiMoon will not run unless the  correct parameters are used. Which parameters are required is based on the dynamics to include flags inside `runprops.txt`.

### Priors, "Object"_priors_df.csv

The priors input file contains a table with the desired priors for your MultiMoon run. As of right now, only uniform priors are well tested and used regularly. This file should have every parameter in it. Thankfully, parameters that appear in the priors file, but aren't being used in the run, will not cause any significant issues. Generally, the prior files are edited only rarely. 

### The Run properties, "runprops.txt"

All of the variables necessary for running MultiMoon are defined in the runprops dictionary. The runprops.txt file is set up in such a way to display the whole dictionary in a readable format, and should look familiar to anyone versed in JSON. A guide for all runprops inputs is located at
