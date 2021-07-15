# MINEOS_synthetics
Calculate mode-based dispersion for layered models using MINEOS:
- phase & group velocity dispersion (Rayleigh and Love)
- sensitivity kernels (Vsv, Vsh, Vpv, Vph, eta, rho; A, C, F, L, N)
- eigenfunctions (toroidal and spheroidal)

Mode tables are used as input to idagrn6, which produces synthetic seismograms for the given layered model. This package also includes branch stripping which produces seismograms for individual modes as well as the excitation terms for each mode.

## Contents
- ./FORTRAN : contains all Fortran binaries required to build MINEOS and idagrn6
- ./run_MINEOS : contains the MATLAB wrappers for running MINEOS to build the mode tables and idagrn6 to calculate synthetic seismograms

## Getting Started

Must have installed
- gfortran (other Fortran compilers might work but have not been tested)
- MATLAB

### For compiling fortran binaries, see [./FORTRAN/README](https://github.com/jbrussell/MINEOS_synthetics/blob/master/FORTRAN/README)

