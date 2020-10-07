# mm_run.py
# 
# The main script that runs MultiMoon
# Darin Ragozzine
# March 27, 2020

"""Run MultiMoon

Inputs:
	filename that is a JSON text file version of run_props dictionary

Outputs:
	diagnostic information about the run

"""
def initializer():
    import os
    import mm_runprops
    os.environ["OMP_NUM_THREADS"] = "1"
    cwd = os.getcwd()
    runs_file = ''
    if 'runs' in cwd:
        runs_file = os.path.basename(os.path.normpath(cwd))
        os.chdir('../../../src')
    import sys
    import numpy as np
    import pandas as pd
    import emcee
    import random
    import h5py
    from tqdm import tqdm
    import mm_init_guess
    import mm_likelihood
    import mm_make_geo_pos
    import mm_priors
    import mm_relast
    import mm_autorun
    import mm_param
    import mm_clustering
    import mm_analysis
    import warnings
    import shutil
    import commentjson as json
    from multiprocessing import Manager
    from csv import writer
    
    #from multiprocessing import Pool, Manager
    from schwimmbad import MultiPool as Pool
    #from mpipool import Pool
    from multiprocessing import Manager
    
    return sys, np, pd, emcee, random, h5py, mm_runprops, mm_init_guess, mm_likelihood, mm_make_geo_pos, mm_priors, mm_relast, mm_autorun, mm_param, mm_clustering, os, mm_analysis, warnings, shutil, json, writer, Manager, tqdm, Pool

# Read in the run props dictionary
# Adds information autogenerated for this specific run
# Checks other information
import warnings
warnings.filterwarnings("ignore", message="divide by zero encountered")
warnings.filterwarnings("ignore", message="Gimbal lock detected")

if __name__ == '__main__':
    sys, np, pd, emcee, random, h5py, mm_runprops, mm_init_guess, mm_likelihood, mm_make_geo_pos, mm_priors, mm_relast, mm_autorun, mm_param, mm_clustering, os, mm_analysis, warnings, shutil, json, writer, Manager, tqdm, Pool = initializer()
    
    pool = Pool()
    
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
        
    manager = Manager()
    best_llhoods = manager.dict()
    best_llhoods['best_llhood'] = -np.inf
    best_llhoods['best_params'] = []

    runprops = mm_runprops.runprops
    runprops['best_llhood'] = -np.inf

    verbose = runprops.get("verbose")
    nwalkers = runprops.get("nwalkers")
    startfromfile = runprops.get("startfromfile")
    nobjects = runprops.get("numobjects")

    name_dict = runprops.get("names_dict")
    objectnames = []
    for i in name_dict.values():
        objectnames.append(i)

# BP TODO: make an option in runprops to start from the end of another run and just append it
# Generate the intial guess for emcee
# starting guess is given by the user as specified in runprops
# and turned into the official initial guess

# LATER TODO: starting -> initial guess function is specificed by user

    guesses = mm_init_guess.mm_init_guess(runprops)	# maybe more args
# ouptut from init_guess is a dataframe with all the desired parameters to be fit
# Getting relevant checking flags from runprops
    dynamicstoincludeflags = runprops.get("dynamicstoincludeflags")
    includesun = runprops.get("includesun")
    com_offset = runprops.get("com_offset")
    paramnames = list(sum(list(guesses), ()))
# Check to make sure that numobjects equals length of dynamics flag
    if len(dynamicstoincludeflags) != runprops.get("numobjects"):
        print("ERROR: Number of objects given in runprops.txt does not match the length of dynamicstoincludeflags")
        sys.exit()

# Now checking each object sequentially
    for i in range(runprops.get("numobjects")):
        if i == 0:
            if dynamicstoincludeflags[i] == "0":
                if not (("mass_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
                if (("sprate_" + str(i+1) in paramnames) or ("j2r2_" + str(i+1) in paramnames) or
            ("spinc_" + str(i+1) in paramnames) or ("splan_" + str(i+1) in paramnames) or
            ("c22r2_" + str(i+1) in paramnames) or ("spaop_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
            elif dynamicstoincludeflags[i] == "1":
                if not (("mass_" + str(i+1) in paramnames) and ("j2r2_" + str(i+1) in paramnames) and
                ("spinc_" + str(i+1) in paramnames) and ("splan_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
                if (("sprate_" + str(i+1) in paramnames) or
            ("c22r2_" + str(i+1) in paramnames) or ("spaop_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
            elif dynamicstoincludeflags[i] == "2":
                if not (("mass_" + str(i+1) in paramnames) and ("j2r2_" + str(i+1) in paramnames) and
                ("spinc_" + str(i+1) in paramnames) and ("splan_" + str(i+1) in paramnames) and
                ("c22r2_" + str(i+1) in paramnames) and ("spaop_" + str(i+1) in paramnames) and
                ("sprate_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
            else:
                print("ERROR: dynamicstoincludeflags contains unallowed numbers. Allowed numbers are 0, 1, 2.")
                sys.exit()
        else:
            if dynamicstoincludeflags[i] == "0":
                if not (("mass_" + str(i+1) in paramnames) and ("sma_" + str(i+1) in paramnames) and
                ("ecc_" + str(i+1) in paramnames) and ("inc_" + str(i+1) in paramnames) and
                ("aop_" + str(i+1) in paramnames) and ("lan_" + str(i+1) in paramnames) and
                ("mea_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
                if (("sprate_" + str(i+1) in paramnames) or ("j2r2_" + str(i+1) in paramnames) or
            ("spinc_" + str(i+1) in paramnames) or ("splan_" + str(i+1) in paramnames) or
            ("c22r2_" + str(i+1) in paramnames) or ("spaop_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
            elif dynamicstoincludeflags[i] == "1":
                if not (("mass_" + str(i+1) in paramnames) and ("sma_" + str(i+1) in paramnames) and
                ("ecc_" + str(i+1) in paramnames) and ("inc_" + str(i+1) in paramnames) and
                ("aop_" + str(i+1) in paramnames) and ("lan_" + str(i+1) in paramnames) and
                ("mea_" + str(i+1) in paramnames) and ("j2r2_" + str(i+1) in paramnames) and
                        ("spinc_" + str(i+1) in paramnames) and ("splan_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
                if (("sprate_" + str(i+1) in paramnames) or
            ("c22r2_" + str(i+1) in paramnames) or ("spaop_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
            elif dynamicstoincludeflags[i] == "2":
                if not (("mass_" + str(i+1) in paramnames) and ("sma_" + str(i+1) in paramnames) and
                ("ecc_" + str(i+1) in paramnames) and ("inc_" + str(i+1) in paramnames) and
                ("aop_" + str(i+1) in paramnames) and ("lan_" + str(i+1) in paramnames) and
                ("mea_" + str(i+1) in paramnames) and ("j2r2_" + str(i+1) in paramnames) and
                        ("spinc_" + str(i+1) in paramnames) and ("splan_" + str(i+1) in paramnames) and
                ("c22r2_" + str(i+1) in paramnames) and ("spaop_" + str(i+1) in paramnames) and
                ("sprate_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
            else:
                print("ERROR: dynamicstoincludeflags contains unallowed numbers. Allowed numbers are 0, 1, 2.")
                sys.exit()
    
    # Now checking the includesun flag
    if includesun:
        if not (("mass_0" in paramnames) and ("sma_0" in paramnames) and
            ("ecc_0" in paramnames) and ("inc_0" in paramnames) and
            ("aop_0" in paramnames) and ("lan_0" in paramnames) and
            ("mea_0" in paramnames)):
            print("ERROR: includesun flag does not match inputs.")
            sys.exit()
    if not includesun:
        if (("mass_0" in paramnames) or ("sma_0" in paramnames) or
        ("ecc_0" in paramnames) or ("inc_0" in paramnames) or
        ("aop_0" in paramnames) or ("lan_0" in paramnames) or
        ("mea_0" in paramnames)):
            print("ERROR: includesun flag does not match inputs.")
            sys.exit()

    # Now checking the com_offset flag
    if com_offset:
        if not (("lat_offset" in paramnames) and ("long_offset" in paramnames)):
            print("ERROR: com_offset flag does not match inputs.")
            sys.exit()
    if not com_offset:
        if (("lat_offset" in paramnames) or ("long_offset" in paramnames)):
            print("ERROR: com_offset flag does not match inputs.")
            sys.exit()

    #ndim is equal to the number of dimension, should this be equal to the number of columns of the init_guess array?
 
    # Convert the guesses into fitting units and place in numpy array
    p0,float_names,fixed_df,total_df_names,fit_scale = mm_param.from_param_df_to_fit_array(guesses,runprops)

    ndim = len(p0[0])
    #we still do not have a constraints or fit scale defined
    
    # Now get observations data frame
    # DS TODO: take observations data frame from runprops
    obsdata = runprops.get('obsdata_file')
    obsdf = 0
    if os.path.exists(obsdata):
        if verbose:
            print("Observational data file " + obsdata + " will be used")
        obsdf = pd.read_csv(obsdata)
    else:
        if verbose:
            print("ERROR: No observational data file exists. Aborting run.")
        sys.exit()

    # Calculating the degrees of freedom
    nobservations = 0
    for i in range(1, nobjects):
        obsdata = obsdf["DeltaLat_" + objectnames[i]].values
        for j in range(len(obsdata)):
            if not np.isnan(obsdata[j]):
                nobservations += 1
        obsdata = obsdf["DeltaLong_" + objectnames[i]].values
        for j in range(len(obsdata)):
            if not np.isnan(obsdata[j]):
                nobservations += 1
    best_llhoods['deg_freedom'] = nobservations - ndim
    
    # Check to see if geocentric_object_position.csv exists and if not creates it
    objname = runprops.get('objectname')
    geofile = '../runs/'+objname+'/'+runprops.get('run_file')+'/geocentric_' + objname + '_position.csv'
    if os.path.exists(geofile):
        if verbose:
            print("File " + geofile + " will be used")
    else:
        if verbose:
            print("No object geocentric position file exists. Creating new file.")
        times = obsdf['time'].tolist()
        mm_make_geo_pos.mm_make_geo_pos(objname, times)	# This is basically a function based on DS's makeHorFile
        if verbose:
            print("geocentric_" + objname + "_position.csv has been created")
    
    # Reads in th geocentric_object data file
    geo_obj_pos = pd.read_csv(geofile)
    
    # Go through initial guesses and check that all walkers have finite posterior probability
    reset = 0
    maxreset = runprops.get("maxreset")
    print('Testing to see if initial params are valid')
    for i in tqdm(range(nwalkers)):  
        llhood = mm_likelihood.log_probability(p0[i,:], float_names,fixed_df.iloc[[i]],total_df_names, fit_scale, runprops, obsdf, geo_obj_pos, best_llhoods)
        reset = 0
        #print(llhood)
        while (llhood == -np.Inf):
            p = random.random()
            p0[i,:] = (p*p0[random.randrange(nwalkers),:] + (1-p)*p0[random.randrange(nwalkers),:])
            llhood = mm_likelihood.log_probability(p0[i,:], float_names,fixed_df,total_df_names, fit_scale, runprops, obsdf,geo_obj_pos, best_llhoods)
            reset += 1
            if reset > maxreset:
                print("ERROR: Maximum number of resets has been reached, aborting run.")
                sys.exit() 

    #import mm_optimize
    #p0 = mm_optimize.mm_optimize(nwalkers, p0, float_names, fixed_df, total_df_names, fit_scale, runprops, obsdf, geo_obj_pos, best_llhoods)
    #sys.exit()

    runprops["is_mcmc"] = True

    # We now have an initial guess for each walker that is not really bad.
    # Begin MCMC
    p0 = list(p0)
    # Now creating the sampler object
    filename = runprops.get("results_folder")+ "/chain.h5"
    
    # BP TODO: make an option in runprops to start from the end of another run and just append it
    
    backend = emcee.backends.HDFBackend(filename)
    backend.reset(nwalkers, ndim)
    moveset = [(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2),]
    moveset = [(emcee.moves.StretchMove(), 1.0),]
    
    the_names = []
    for i in total_df_names:
        the_names.append(i[0])
    
    if runprops.get('updatebestfitfile'):
        the_file = runprops.get('results_folder') + '/best_likelihoods.csv'
        with open(the_file, 'a+', newline='') as write_obj:
            csv_writer = writer(write_obj, delimiter = ',')
            the_names.insert(0,'Prior')
            the_names.insert(0,'Reduced chi-sq')
            the_names.insert(0,'Likelihood')
            for i in range(runprops.get('numobjects')-1):
                the_names.append('Residuals_Lon_Obj_'+str(i+1))
                the_names.append('Residuals_Lat_Obj_'+str(i+1))
            csv_writer.writerow(the_names)
            
            
        
    with Pool(runprops.get("numprocesses")) as pool:
    #with Pool() as pool:
        
        if not pool.is_master():
            pool.wait()
            sys.exit(0)
    
        sampler = emcee.EnsembleSampler(nwalkers, ndim, 
        mm_likelihood.log_probability, backend=backend, pool = pool,
            args = (float_names, fixed_df, total_df_names, fit_scale, runprops, obsdf,geo_obj_pos, best_llhoods),
            moves = moveset)
        print('sampler created')

    #Starting the burnin
    # BP TODO: autoburnin??
    # So looking at how the emcee documentation does burn ins while saving the file, it seems like
    # the best way to do a run is to just do a single long run until the ess > 100 and cut the 
    # burn in off afterwards. This way you can save the whole chain and you can really analyze where to
    # cut off the burn in.
    # I think i want to still create an autoburnin but I really would like to look at a completed
    # run to see what the burn in looks like... It should be a few autocorrelation times
    
        nburnin = runprops.get("nburnin")
        if verbose:
            print("Starting the burn in")
    
        state = sampler.run_mcmc(p0, nburnin, progress = True, store = True)

        # Now running the clustering algorithm! (if desired)
        if runprops.get("use_clustering") and nburnin != 0:
            sampler, state = mm_clustering.mm_clustering(sampler, state, float_names, fixed_df, total_df_names, fit_scale, runprops, obsdf,geo_obj_pos, best_llhoods, backend, pool, mm_likelihood, ndim, moveset)
        
        sampler.reset()

    # Now do the full run with essgoal and initial n steps
    
        nsteps = runprops.get("nsteps")
        essgoal = runprops.get("essgoal")
        maxiter = runprops.get("maxiter")
        initsteps = runprops.get("nsteps")
        
        sampler,ess = mm_autorun.mm_autorun(sampler, essgoal, state, initsteps, maxiter, verbose, objname, p0, runprops)
        
    print("effective sample size = ", ess)
    chain = sampler.get_chain(thin = runprops.get("nthinning"))
    flatchain = sampler.get_chain(flat = True, thin = runprops.get("nthinning"))
        
    # Begin analysis!
    print('Beginning mm_analysis plots')
    
    mm_analysis.plots(sampler, guesses.columns, objname, fit_scale, float_names, obsdf, runprops, geo_obj_pos, mm_make_geo_pos)
    runpath = runprops.get("results_folder")+"/runprops.txt"
    
    with open(runpath, 'w') as file:
        file.write(json.dumps(runprops, indent = 4))
        
    if runprops.get('build_init_from_llhood'):
        csvfile = runprops.get("resuts_folder")+"/best_likelihoods.csv"
        likelihoods = pd.read_csv(csvfile, sep = '\t', header = 0)
        
        
        params = likelihoods.iloc[[-1]].transpose()
        params = params.drop(['Likelihood'], axis=0)
        for i in range(runprops.get('numobjects')-1):
            params = params.drop(['Residuals_Lat_Obj_'+str(i+1),'Residuals_Lon_Obj_'+str(i+1)], axis =0)

        init_guess = pd.read_csv(runprops.get('init_filename'))
        stddev  = init_guess['stddev'].tolist()

        params['stddev'] = stddev
        params.columns = ['mean', 'stddev']
        #print(params)
        new_init = "../data/"+runprops.get("objectname")+"/"+runprops.get("objectname")+"_init_guess_from_llhood.csv"
        params.to_csv(new_init, sep = ',')
        

    # make other diagnostic plots
    # TODO: orbit astrometry plots
    # TODO: residual plots
    
def run():
    sys, np, pd, emcee, random, h5py, mm_runprops, mm_init_guess, mm_likelihood, mm_make_geo_pos, mm_priors, mm_relast, mm_autorun, mm_param, mm_clustering, os, mm_analysis, Pool, warnings, shutil, json, writer, Manager = initializer()
    
    manager = Manager()
    best_llhoods = manager.dict()
    best_llhoods['best_llhood'] = -np.inf
    best_llhoods['best_params'] = []

    runprops = mm_runprops.runprops
    runprops['best_llhood'] = -np.inf

    verbose = runprops.get("verbose")
    nwalkers = runprops.get("nwalkers")
    startfromfile = runprops.get("startfromfile")
    nobjects = runprops.get("numobjects")

    name_dict = runprops.get("names_dict")
    objectnames = []
    for i in name_dict.values():
        objectnames.append(i)

# BP TODO: make an option in runprops to start from the end of another run and just append it
# Generate the intial guess for emcee
# starting guess is given by the user as specified in runprops
# and turned into the official initial guess

# LATER TODO: starting -> initial guess function is specificed by user

    guesses = mm_init_guess.mm_init_guess(runprops)	# maybe more args
# ouptut from init_guess is a dataframe with all the desired parameters to be fit

# Getting relevant checking flags from runprops
    dynamicstoincludeflags = runprops.get("dynamicstoincludeflags")
    includesun = runprops.get("includesun")
    paramnames = list(sum(list(guesses), ()))
# Check to make sure that numobjects equals length of dynamics flag
    if len(dynamicstoincludeflags) != runprops.get("numobjects"):
        print("ERROR: Number of objects given in runprops.txt does not match the length of dynamicstoincludeflags")
        sys.exit()

# Now checking each object sequentially
    for i in range(runprops.get("numobjects")):
        if i == 0:
            if dynamicstoincludeflags[i] == "0":
                if not (("mass_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
                if (("sprate_" + str(i+1) in paramnames) or ("j2r2_" + str(i+1) in paramnames) or
            ("spinc_" + str(i+1) in paramnames) or ("splan_" + str(i+1) in paramnames) or
            ("c22r2_" + str(i+1) in paramnames) or ("spaop_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
            elif dynamicstoincludeflags[i] == "1":
                if not (("mass_" + str(i+1) in paramnames) and ("j2r2_" + str(i+1) in paramnames) and
                ("spinc_" + str(i+1) in paramnames) and ("splan_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
                if (("sprate_" + str(i+1) in paramnames) or
            ("c22r2_" + str(i+1) in paramnames) or ("spaop_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
            elif dynamicstoincludeflags[i] == "2":
                if not (("mass_" + str(i+1) in paramnames) and ("j2r2_" + str(i+1) in paramnames) and
                ("spinc_" + str(i+1) in paramnames) and ("splan_" + str(i+1) in paramnames) and
                ("c22r2_" + str(i+1) in paramnames) and ("spaop_" + str(i+1) in paramnames) and
                ("sprate_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
            else:
                print("ERROR: dynamicstoincludeflags contains unallowed numbers. Allowed numbers are 0, 1, 2.")
                sys.exit()
        else:
            if dynamicstoincludeflags[i] == "0":
                if not (("mass_" + str(i+1) in paramnames) and ("sma_" + str(i+1) in paramnames) and
                ("ecc_" + str(i+1) in paramnames) and ("inc_" + str(i+1) in paramnames) and
                ("aop_" + str(i+1) in paramnames) and ("lan_" + str(i+1) in paramnames) and
                ("mea_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
                if (("sprate_" + str(i+1) in paramnames) or ("j2r2_" + str(i+1) in paramnames) or
            ("spinc_" + str(i+1) in paramnames) or ("splan_" + str(i+1) in paramnames) or
            ("c22r2_" + str(i+1) in paramnames) or ("spaop_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
            elif dynamicstoincludeflags[i] == "1":
                if not (("mass_" + str(i+1) in paramnames) and ("sma_" + str(i+1) in paramnames) and
                ("ecc_" + str(i+1) in paramnames) and ("inc_" + str(i+1) in paramnames) and
                ("aop_" + str(i+1) in paramnames) and ("lan_" + str(i+1) in paramnames) and
                ("mea_" + str(i+1) in paramnames) and ("j2r2_" + str(i+1) in paramnames) and
                        ("spinc_" + str(i+1) in paramnames) and ("splan_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
                if (("sprate_" + str(i+1) in paramnames) or
            ("c22r2_" + str(i+1) in paramnames) or ("spaop_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
            elif dynamicstoincludeflags[i] == "2":
                if not (("mass_" + str(i+1) in paramnames) and ("sma_" + str(i+1) in paramnames) and
                ("ecc_" + str(i+1) in paramnames) and ("inc_" + str(i+1) in paramnames) and
                ("aop_" + str(i+1) in paramnames) and ("lan_" + str(i+1) in paramnames) and
                ("mea_" + str(i+1) in paramnames) and ("j2r2_" + str(i+1) in paramnames) and
                        ("spinc_" + str(i+1) in paramnames) and ("splan_" + str(i+1) in paramnames) and
                ("c22r2_" + str(i+1) in paramnames) and ("spaop_" + str(i+1) in paramnames) and
                ("sprate_" + str(i+1) in paramnames)):
                    print("ERROR: dynamics to include flags does not match the input parameters for object " + str(i+1))
                    sys.exit()
            else:
                print("ERROR: dynamicstoincludeflags contains unallowed numbers. Allowed numbers are 0, 1, 2.")
                sys.exit()
    
    # Now checking the includesun flag
    if includesun:
        if not (("mass_0" in paramnames) and ("sma_0" in paramnames) and
            ("ecc_0" in paramnames) and ("inc_0" in paramnames) and
            ("aop_0" in paramnames) and ("lan_0" in paramnames) and
            ("mea_0" in paramnames)):
            print("ERROR: includesun flag does not match inputs.")
            sys.exit()
    if not includesun:
        if (("mass_0" in paramnames) or ("sma_0" in paramnames) or
        ("ecc_0" in paramnames) or ("inc_0" in paramnames) or
        ("aop_0" in paramnames) or ("lan_0" in paramnames) or
        ("mea_0" in paramnames)):
            print("ERROR: includesun flag does not match inputs.")
            sys.exit()
        
    #ndim is equal to the number of dimension, should this be equal to the number of columns of the init_guess array?
 
    # Convert the guesses into fitting units and place in numpy array
    p0,float_names,fixed_df,total_df_names,fit_scale = mm_param.from_param_df_to_fit_array(guesses,runprops)
    
    
    ndim = len(p0[0])
    #we still do not have a constraints or fit scale defined
    
    # Now get observations data frame
    # DS TODO: take observations data frame from runprops
    obsdata = runprops.get('obsdata_file')
    
    obsdf = 0
    if os.path.exists(obsdata):
        if verbose:
            print("Observational data file " + obsdata + " will be used")
        obsdf = pd.read_csv(obsdata)
    else:
        if verbose:
            print("ERROR: No observational data file exists. Aborting run.")
        sys.exit()

    # Calculating the degrees of freedom
    nobservations = 0
    for i in range(1, nobjects):
        obsdata = obsdf["DeltaLat_" + objectnames[i]].values
        for j in range(len(obsdata)):
            if not np.isnan(obsdata[j]):
                nobservations += 1
        obsdata = obsdf["DeltaLong_" + objectnames[i]].values
        for j in range(len(obsdata)):
            if not np.isnan(obsdata[j]):
                nobservations += 1
    best_llhoods['deg_freedom'] = nobservations - ndim
    
    # Check to see if geocentric_object_position.csv exists and if not creates it
    objname = runprops.get('objectname')
    if os.path.exists("../data/" + objname + "/geocentric_" + objname + "_position.csv"):
        if verbose:
            print("Object geocentric position file geocentric_" + objname + "_position.csv will be used")
    else:
        if verbose:
            print("No object geocentric position file exists. Creating new file.")
        times = obsdf['time'].tolist()
        mm_make_geo_pos.mm_make_geo_pos(objname, times)	# This is basically a function based on DS's makeHorFile
        if verbose:
            print("geocentric_" + objname + "_position.csv has been created")
    
    # Reads in th geocentric_object data file
    geo_obj_pos = pd.read_csv("../data/" + objname + "/geocentric_" + objname + "_position.csv")
    
    # Go through initial guesses and check that all walkers have finite posterior probability
    reset = 0
    maxreset = runprops.get("maxreset")
    
    print('Testing to see if initial params are valid')
    for i in range(nwalkers):  
        llhood = mm_likelihood.log_probability(p0[i,:], float_names,fixed_df.iloc[[i]],total_df_names, fit_scale, runprops, obsdf, geo_obj_pos, best_llhoods)
        reset = 0
        #print(llhood)
        while (llhood == -np.Inf):
            p = random.random()
            p0[i,:] = (p*p0[random.randrange(nwalkers),:] + (1-p)*p0[random.randrange(nwalkers),:])
            llhood = mm_likelihood.log_probability(p0[i,:], float_names,fixed_df,total_df_names, fit_scale, runprops, obsdf,geo_obj_pos, best_llhoods)
            reset += 1
            print(llhood)
            if reset > maxreset:
                print("ERROR: Maximum number of resets has been reached, aborting run.")
                sys.exit() 
    runprops["is_mcmc"] = True

    # We now have an initial guess for each walker that is not really bad.
    # Begin MCMC
    p0 = list(p0)
    # Now creating the sampler object
    filename = "../runs/" + runprops.get("objectname") + "_" + runprops.get("date") + "/chain.h5"
    
    # BP TODO: make an option in runprops to start from the end of another run and just append it
    
    backend = emcee.backends.HDFBackend(filename)
    backend.reset(nwalkers, ndim)
    moveset = [(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2),]
    moveset = [(emcee.moves.StretchMove(), 1.0),]
    
    the_names = []
    for i in total_df_names:
        the_names.append(i[0])
    
    if runprops.get('updatebestfitfile'):
        the_file = runprops.get('runs_folder') + '/best_likelihoods.csv'
        with open(the_file, 'a+', newline='') as write_obj:
            csv_writer = writer(write_obj, delimiter = ',')
            the_names.insert(0,'Prior')
            the_names.insert(0,'Reduced chi-sq')
            the_names.insert(0,'Likelihood')
            for i in range(runprops.get('numobjects')-1):
                the_names.append('Residuals_Lon_Obj_'+str(i+1))
                the_names.append('Residuals_Lat_Obj_'+str(i+1))
            csv_writer.writerow(the_names)
            
        
    with Pool(runprops.get("numprocesses")) as pool:
    
        sampler = emcee.EnsembleSampler(nwalkers, ndim, 
        mm_likelihood.log_probability, backend=backend, pool = pool,
            args = (float_names, fixed_df, total_df_names, fit_scale, runprops, obsdf,geo_obj_pos, best_llhoods),
            moves = moveset)
        print('sampler created')

    #Starting the burnin
    # BP TODO: autoburnin??
    # So looking at how the emcee documentation does burn ins while saving the file, it seems like
    # the best way to do a run is to just do a single long run until the ess > 100 and cut the 
    # burn in off afterwards. This way you can save the whole chain and you can really analyze where to
    # cut off the burn in.
    # I think i want to still create an autoburnin but I really would like to look at a completed
    # run to see what the burn in looks like... It should be a few autocorrelation times
    
        nburnin = runprops.get("nburnin")
        if verbose:
            print("Starting the burn in")
    
        state = sampler.run_mcmc(p0, nburnin, progress = True, store = True)

        # Now running the clustering algorithm! (if desired)
        print(runprops.get("nburnin"))
        if runprops.get("use_clustering") and runprops.get("nburnin") != 0:
            sampler, state = mm_clustering.mm_clustering(sampler, state, float_names, fixed_df, total_df_names, fit_scale, runprops, obsdf,geo_obj_pos, best_llhoods, backend, pool, mm_likelihood, ndim, moveset)
        
        sampler.reset()

    # Now do the full run with essgoal and initial n steps
    
        nsteps = runprops.get("nsteps")
        essgoal = runprops.get("essgoal")
        maxiter = runprops.get("maxiter")
        initsteps = runprops.get("nsteps")
        
        sampler,ess = mm_autorun.mm_autorun(sampler, essgoal, state, initsteps, maxiter, verbose, objname, p0, runprops)
        
    print("effective sample size = ", ess)
    chain = sampler.get_chain(thin = runprops.get("nthinning"))
    flatchain = sampler.get_chain(flat = True, thin = runprops.get("nthinning"))
        
    # Begin analysis!
    print('Beginning mm_analysis plots')
    
    mm_analysis.plots(sampler, guesses.columns, objname, fit_scale, float_names, obsdf, runprops, geo_obj_pos, mm_make_geo_pos)
    runpath = "../runs/"+runprops.get("objectname")+"_"+runprops.get("date")+"/runprops.txt"
    
    with open(runpath, 'w') as file:
        file.write(json.dumps(runprops, indent = 4))
        
    if runprops.get('build_init_from_llhood'):
        csvfile = "../runs/"+runprops.get("objectname")+"_"+runprops.get("date")+"/best_likelihoods.csv"
        likelihoods = pd.read_csv(csvfile, sep = '\t', header = 0)
        
        
        params = likelihoods.iloc[[-1]].transpose()
        params = params.drop(['Likelihood'], axis=0)
        for i in range(runprops.get('numobjects')-1):
            params = params.drop(['Residuals_Lat_Obj_'+str(i+1),'Residuals_Lon_Obj_'+str(i+1)], axis =0)

        init_guess = pd.read_csv(runprops.get('init_filename'))
        stddev  = init_guess['stddev'].tolist()

        params['stddev'] = stddev
        params.columns = ['mean', 'stddev']
        #print(params)
        new_init = "../data/"+runprops.get("objectname")+"/"+runprops.get("objectname")+"_init_guess_from_llhood.csv"
        params.to_csv(new_init, sep = ',')
