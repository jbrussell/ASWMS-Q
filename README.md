# Automated Surface-Wave Measurement System (ASWMS) + anisotropy + attenuation

##### This project builds from the original ASWMS package of Ge Jin ([https://github.com/jinwar/matgsdf](https://github.com/jinwar/matgsdf)) but contains some significant modifications/additions.

From the original repository:
>Ge Jin, James Gaherty
>
>Lamont Doherty Earth Observatory, Columbia University
>
>This program can automatically measure the surface-wave phase velocity of a certain region based on the cross-correlation measurements between the nearby stations and Helmholtz tomography. 
>
>Please read the manual for more details.

The primary add-ons in this package are the ability to solve explicitly for **attenuation** and a new method for solving for **azimuthal anisotropy**.

It is recommended that the waveform download step is performed using this python package: [https://github.com/jbrussell/fetch_EVENTS](https://github.com/jbrussell/fetch_EVENTS) and data are loaded externally via **a1_a_sac2eventmat.m**. While ASWMS has a built in data download and preprocessing script, it does not pull the moment tensor information. We have implemented some new quality control efforts, which now require information about earthquake magnitude (and in some cases, focal mechanism). Most of the scripts should run with the old data download step, but your mileage may vary.

### Attenuation

Rayleigh-wave attenuation and site amplification are solved for following Bao et al. (2016) GJI [doi: 10.1093/gji/ggw151](https://academic.oup.com/gji/article/206/1/241/2606522) and Lin et al. (2012) JGR [doi: doi:10.1029/2012JB009208](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2012JB009208) in **b1_estimate_alpha_beta_Bao16_errs_bs_gammacorr.m**. This requires both travel-time and amplitude fields and their spatial derivatives. We have included the option to apply station corrections following Eddy and Ekström [doi: 10.1016/j.epsl.2014.01.013](https://www.sciencedirect.com/science/article/pii/S0012821X14000223?via%3Dihub) in **a8_a_receiver_terms.m**, under the assumption of an even azimuthal distribution of events, but these station corrections **SHOULD NOT BE APPLIED** and should only be used to compare to site amplification, β.

The calculation requires gradients and Laplacian of the travel-time and amplitude fields. We've included two ways to estimate those fields. The first follows the original methods of fitting a surface to the amplitudes or traveltimes and simply taking the first and second spatial derivatives. However, we have found that this is not always a stable calculation, particularly when station geometry is sparse. Additionally, it is difficult to control the smoothness of the final Laplacian map by smoothing the original surface. 

The second (new) method utilizes the ray tomography inversion implemented by Jin & Gaherty (2015) GJI [doi: 10.1093/gji/ggv079](https://doi.org/10.1093/gji/ggv079) which uses interstation travel times to solve for the spatial gradient field of travel time. We apply an identical inversion approach to the amplitude data in order to solve also for the amplitude gradient field (**a8_ampgrad_inv.m** and **a8_ampgrad_norm_inv.m**) and simply take the first spatial derivatives of those maps to get the Laplacian. The advantage of this approach is that the regularized inversion provides a more stable measurement. In addition, the smoothness of the Laplacian field can be directly controlled by the 2nd derivative smoothing constraint in the inversion.

### Azimuthal Anisotropy (*WORK IN PROGRESS*)

Azimuthal anisotropy can be solved for using two different techniques. The first technique, which was implemented in the original ASWMS, package attempts to solve for anisotropy in the eikonal stacking phase, after the phase velocity maps have already been estimated (**a7_b_stack_phv_aniso_1D.m** or **a7_b_stack_phv_aniso_2D.m**). Azimuthal anisotropy is inferred at each grid point by gathering observed eikonal phase velocities from neighboring grid points and fitting a 2-theta sinusoid to residuals with respect to the mean. We have also extended this to the helmholtz stacking phase (**a9_c_stack_helm_aniso_1D.m** or **a9_c_stack_helm_aniso_2D.m**).

The second (new) method is to invert for azimuthal anisotropy and phase velocity simultaneously using all earthquakes at once (**a6_b_eikonal_eq_2DanisoRT.m** or **a6_b_eikonal_eq_flat_1DanisoRT.m**). The advantage is that both isotropic velocity and anisotropy are parameterized explicitly in the inversion. This single-step solution utilizes all data at once so outliers and bad measurements should be less problematic, in theory. In addition, for the single-step the smoothing you apply is the smoothing you get in the final map, whereas for the per-event inversion you control only smoothing of each individual event map but the final stacked map can end up much smoother (in my experience). The main disadvantage is that the Helmholtz correction cannot be applied within this framework. 

There are options for both a 1-D anisotropy (one anisotropy strength and fast azimuth for the entire array) and 2-D anisotropy (anisotropy strength and azimuth at each grid point).

___
Please cite:

Russell, J. B., & Dalton, C. A. (2022). Rayleigh wave attenuation and amplification measured at ocean-bottom seismometer arrays using Helmholtz tomography. Journal of Geophysical Research: Solid Earth, 127, e2022JB025174. [https://doi.org/10.1029/2022JB025174](https://doi.org/10.1029/2022JB025174)

Jin, G., & Gaherty, J. B. (2015). Surface wave phase-velocity tomography based on multichannel cross-correlation. Geophysical Journal Interna-
tional, 201(3), 1383–1398. [https://doi.org/10.1093/gji/ggv079](https://doi.org/10.1093/gji/ggv079)
