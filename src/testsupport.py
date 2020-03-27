# testsupport.py
#
# Darin Ragozzine
# March 27, 2020
# 
# Supports testing of MultiMoon by including functions
# that return valid "shapes" and types of various
# planned pieces of MultiMoon. 


# test_runprops
# 

def test_runprops():
    """ Returns a runprops dictionary. """
    
    runprops = {
        'date': 'today',
        'time': 'now',
        'user': 'autogenerate'
        'MultiMoon commit hash': 'autogenerate'
        'mm_initguess_which' : 'get_init_guess_from_dist' # which function from mm_initguess to use
        'numobjects' : 3
        'RunGoal' 'user puts their goal for this run here'
        'objectname' : 'Haumea'
        'nwalkers' : 30
        'nsteps' : 1001
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
