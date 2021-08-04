Josh Russell - July 20, 2021
This is a new version of MINEOS from CIG [https://geodynamics.org/cig/software/mineos/]. Colleen Dalton modified this version slightly to output .eig files that are compatible with all our other MINEOS fortran codes. I further modified the inputs slightly in order to be perfectly compatible with our Matlab wrapper.


Note from Colleen:
"I’m emailing you since I know you are regular users of mineos and other programs that utilize the output of mineos. I have spent the last few days playing with the CIG version of mineos. It seems more stable than the legacy version many (all?) of us have been using. Notably, when running a full mode catalog, it does not get hung up so easily. The CIG distribution generates .asc. and .eig files, just like our version does, but the binary .eig file is formatted differently and therefore not compatible with mineos_q, mineos_strip, mineos_table, and so on. 

I therefore modified the CIG version so that the .eig file now is compatible with those other programs (see attached). For once I carefully documented all changes with my initials (CAD). I figured I’d share it with you in case you find it helpful. Here are a few pros and cons of the (CAD-modified) CIG version: 

Pros:
- seems to not get hung up easily on overtones
- all routines are contained with a single .f file that can be compiled with gfortran; also, no parameter.h file
- for a single Earth model (ATL2a), I have tested the spheroidal and toroidal output against our legacy version. I checked that the phase velocity, group velocity, Q, Frechet kernels, and synthetic seismograms all agree (they do!)
- the .asc file contains phase velocity so no need to run mineos_q just for that

Cons: 
- the variable for max number of rows (mk) is hard-wired in numerous routines, so if you want to change it, you need to make sure you change it in all subroutines (easily done with find and replace) 
- I have not tested additional Earth models
- I have not tested radial modes
- There are a couple of changes I needed to make that I’m not sure why they were necessary. For the spheroidal eigenfunctions,  I needed to multiply their values by -1. And the Q values calculated for the toroidal modes do not seem correct; I found I needed to run mineos_q on the toroidal modes to get correct Q values."