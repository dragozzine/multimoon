# Multimoon–SPINNY README
*****************************************************************
Author: Seth Pincock

Contributors: Darin Ragozzine, Brigham Young University

Acknowledgements: Simon Porter, Southwest Research Institute
*****************************************************************
  This program is coded using Python 3.5.2. SPINNY (the integrator itself) is coded in C++ and Cythonized.
  
  The SPINNY (SPIN+N-bodY) integrator was developed by Dr. Simon Porter, and is designed to work for any N-body system, with the Sun included. This program was written by Seth Pincock in order to efficiently utilize the SPINNY integrator as a part of the collaborative multimoon project, with the goal of producing accurate fits for the non-Keplerian orbits of N-body Kuiper Belt Object systems, including the quadrupole shapes of the bodies. 
  
The program has a number of automated functions:
1. Perform an N-body integration for a system over an inputted list of times
    * If a keplerian, 2-body system is detected, a simpler 2-body integration is run instead for efficiency.
2. Generate plots of the output data from a SPINNY integration  
3. Generate a 3D-rendered animation of the system, using VPython
    * An animation of either the bodies' orbits or their spins can be generated
  
**NOTE:** This program is not yet capable of computing accurately the spin dynamics of the bodies due to unresolved issues within SPINNY. As such, it is only capable of accurately computing point masses. The spin data are still included in the integrations, but the output spin data are not physically accurate (as of May, 2020).
*****************************************************************
## Inputs/Outputs
**ALL input data MUST be in units of kilometers, kilograms, seconds, and radians**

Input:
* Physical parameters for each of the bodies 
  * Mass
  * J<sub>2</sub>R<sup>2</sup>
  * C<sub>22</sub>R<sup>2</sup>
  * Length of the smallest semi-axis (c)
* Orbital parameters for each of the bodies 
  * Semi-major axis
  * Eccentricity
  * Incination
  * Longitude of the ascending node
  * Argument of pericenter
  * Mean anomaly at epoch
* Spin Parameters (NOTE: not currently working, but support still exists in the program's input)
  * Angular velocity
  * Precession angle
  * Obliquity angle
  * Longitudinal angle

*If no input is given for any of these values, a default value is automatically assigned.*  
    
Output:
* x,y,z POSITION components for each body
* x,y,z VELOCITY components for each body
(position and velocities are given with respect to the primary)
* Orbital parameters
  * Inclination (with respect to the ecliptic)
  * Eccentricity
  * Longitude of the ascending node
  * Argument of periapsis
  * Mean anomaly
  * Semi-major axis 
  * (Future versions will include spin parameters and quaternions)
* Spin Parameters (**NOTE: Functionality of spins still incomplete**)
  * Angular velocity
  * Precession angle
  * Obliquity angle


The spin orientation angles (precession, obliquity, and longitude) correspond to the longitude of the ascending node, inclination, and argument of pericenter, respectively. The angles also correspond to the ZXZ Euler transformation. **The orientation angles should be defined in the inertial (ecliptic) reference frame, NOT the orbit frame.**

****************************************************************
## Prerequisites

In order to run, the program requires the following files:
- SPINNY files:
  - cashkarp.hpp
  - conics.hpp
  - spinny.hpp
  - spinny.pyx
  - setup.py
- Multimoon files:
  - keplerian.py
  - mm_vpython.py
  - ~~quaternion.py~~  (no longer vital? These functions replaced with SciPy spatial transform)
  - spinny_generate.py
  - spinny_plots.py
  - **SPINNY_User_Interface.py** (the only script you should need to execute directly)

****************************************************************
## Using SPINNY:

1. If it has not been already, SPINNY will need to be installed. SPINNY is a C++ code that will need to be Cythonized in order to be read by the python scripts. To do this, execute one of the following from a command line:
    1. `python setup.py build_ext --inplace` to install in the local directory (useful for debugging. This command will need to be re-executed every time spinny.hpp or spinny.pyx is edited.)
    2. `python setup.py install --user` to install for the current user
    3. `sudo python setup.py install` to install globally
2. Once SPINNY is installed, a new, computer-generated file called "spinny.cpp" should appear in the local directory. If SPINNY was already installed, the old "spinny.cpp" will be overwritten when "setup.py" is re-run. 
3. From a command line, execute `python SPINNY_User_Interface.py` or from a Jupyter notebook `run SPINNY_User_Interface`. 
4. The program is highly automated, and by design requires little user input. Follow the instructions printed onscreen in order to begin the desired function. **Note that any input files MUST be .csv files.**
5. Any output files will be saved to the local directory. Note that the plots may not be visible if run from within a terminal window, but will still save to the directory. If run from a Jupyter notebook, run `%matplotlib inline` before the Interface is run in order to ensure the figures will be visible within the notebook environment.

  


