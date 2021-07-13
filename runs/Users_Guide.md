# A User's Guide to Using MultiMoon

### Purpose of MultiMoon

MultiMoon is a program designed to fit solar system small body multiple systems with non-Keplerian orbits including their quadrupole shapes. This will allow us to learn more about distant solar system objects, and determine their orbital parameters through Bayesian statistical sampling. 
This code requires several types of files to be able to run, which will be exaplined here. In particular, to get MultiMoon to functino properly, you will need to provide a...
1) A csv file of the Observational astrometry of the object
2) A csv file of the inital guess for the physical and orbital parameters of the object
3) A csv file that will control how the Priors of the given parameters will be calculated.

Additonally a .txt file called 'runprops.txt' will need to be initialized which holds all of the variables needed to control the input and output for the program. This file will be explained in depth, and the variables explained. 

In the exmaples for csv file organiation, the Haumea system (Primary = Haumea, Satellites = Hiiaka, Namaka) will be the example for the set up. 



### Observation CSV,  "Object"_obs_df.csv

This file represents a dataframe that is read into MultiMoon and then used as the posterior belief against which our model will be calculated. As a dataframe, it will have named columns. As a dataframe the columns do not need to be in order, besides the first column, which will represent the indices, which in this case is 'time', and the indice will portray the individual epoch's of the observations listed. The column names often include the name of the object and the satellites. The columns should be organized like so....

###### time, Lat_Prim, Long_Prim, DeltaLat_Hiiaka, DeltaLong_Hiiaka, DeltaLat_Hiiaka_err, DeltaLong_Hiiaka_err, DeltaLat_Namaka, DeltaLong_Namaka, DeltaLat_Namaka_err, DeltaLong_Namaka_err

The astrometry in this file should all be in the J2000 ecliptic. If current data is not given in the J2000 ecliptic, but rather the equatorial place, a function is currently being written for conversion to the ecliptic plane. The order of the columns does not matter, since the columns are called through the use of pandas dataframes, which uses the names of the satellites. The names must be exact and must match the names_dictionary in the runprops dictionary, otherwise the code will throw an error.

Observations for objects at certain times can be left as empty, if, for example, the satellite is located behind the primary at the time of observations. This will be translated by the code as a NaN observation at that time, and the code will compensate for that non-observation.

The observation dataframe that you would like to use must be specified in the runprops.txt file as the obsdata_file variable.

### Initial Guess, "Object"_init_guess.csv

The initial guess file is composed of two columns, and a number of rows equal to the numebr of parameters that the user wants to test, which each row representing a parameter for which we are initializing our code. The dynamicstoincludeflags variable in the runrops dictinoary controls which variables must be included in the dataframe, and will be described later. 

The first column asks for the mean of each parameter. In most cases, this will be the currently best known value for that parameter. The second column will ask for the standard deviation of error that you would like to initalize the code with. The inital guess csv will be used to create a number of variables for each floating parameter and will use a normal distribution to spread out the initial sample of parameters.  


### Priors Dataframe, "Object"_priors_df.csv

The priors csv file is meant to control how we calculate our prior beliefs for our Bayesian statistical calculations. Each row, like the inital guess cv file, represents a parameter, floating or fixed. The priors dataframe does not require you to include a row for every parameter, unlike the initial guess dataframe. 

The first column of the file is called dist_shape. This column should include an integer from 0-3. Any other number will cause this parameter's prior belief to be skipped during the calculation of the total prior likelihood. Each integer 0-3 means a different distribution shape will be used to calculate that parameter's prior likelihood.

0 - A Uniform Distirbution will be calculated
1 - A Log-Uniform Distribution will be calculated
2 - A Normal Distribution will be calculated
3 - A Log-Normal Distribution will be calculated

The remaining columns control the distribution calculations. 
Columns 2 and 3 are uni-low and uni-up, representing the lower and upper limits of the desired uniform distirbution calulcation.
Columns 4 and 5 are log-uni-low and log-uni-up, representing the lower and upper limits of the desired uniform distirbution calulcation.
Columns 6 and 7 are norm-cen and norm-spread, representing the center of the normal distiribution and the standard deviation expected for the parameter
Columns 8 and 9 are log-norm-cen and log-norm-spread, representing the center of the log-normal distiribution and the standard deviation expected for the parameter


### The Run properties, "runprops.txt"

All of the variables necessary for running MultiMoon are defined in the runprops dictionary. The runprops.txt file is set up in such a way to display the whole dictionary in a readable format, and should look familiar to anyone versed in JSON.

Some variables are very necessary for MultiMoon to run, and others are just cosmetic, and help keep a record of the goals and pruposes of the run, for future reference. 
