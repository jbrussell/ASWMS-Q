# Two-plane wave
###### jbrussell - 10/18/2021

This directory contains a set of codes for converting ASWMS measurements to the format required for the two-plane wave (TPW) inversion scheme of Forsyth & Li (2005); Yang & Forsyth (2006). The main pieces of the inversion are written in Fortran, and here we've written a more user friendly Matlab wrapper around those scripts. 

There are many different version of the TPW that have been developed over the years. This version is based on srchwave589, which was later modified slightly by Zhitu Ma for the NoMelt project. The most complete (up-to-date?) list of the different versions of TPW and their differences are here [https://github.com/chukren/SurfwaveTomoPrograms/tree/master/SurfaceProgs/phaseVelInversion](https://github.com/chukren/SurfwaveTomoPrograms/tree/master/SurfaceProgs/phaseVelInversion).

This version of TPW inverts for: 
- 2D phase velocity
- 1D azimuthal anisotropy (can solve for 2D anisotropy by setting `iarea = nnodes` in `./fortran/srchwave589.JdF.nophase2.iarea.jbr.f`, I think)
- 1D attenuation parameter (alpha)
- Station amplitude correction terms

### Setting up directories

After running ASWMS to at least `a4_gsdfmain.m`, copy the desired ASWMS working directory at location `parameters.workingdir` to the `./TPW` directory. Then, open `setup_parameters_tpw.m` and modify `path2ASWMS_output` to point to this new working directory location.

### Compiling the Fortran code
```
cd ./fortran
make
cd ..
```