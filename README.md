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

The primary add-ons in this package are the ability to solve explicitly for **azimuthal anisotropy** and **attenuation**.

It is recommended that the waveform download step is performed using this python package: [https://github.com/jbrussell/fetch_EVENTS](https://github.com/jbrussell/fetch_EVENTS) and data are loaded externally via **a1_a_sac2eventmat.m**. While ASWMS has a built in data download and preprocessing script, it does not pull the moment tensor information including earthquake magnitude. We have implemented some new quality control efforts, which now require this information. Most of the scripts should run with the old data download step, but your mileage may vary.

### Azimuthal Anisotropy (*WORK IN PROGRESS*)

Azimuthal anisotropy can be solved for using two different techniques. The first technique which was implemented in the original ASWMS package attempts to solve for anisotropy in the eikonal stacking phase, after the phase velocity maps have already been estimated (**a7_b_stack_phv_aniso.m**). Azimuthal anisotropy is inferred at each grid point by gathering observed eikonal phase velocities from neighboring grid points and fitting a 2-theta sinusoid to residuals with respect to the mean.

The second (new) method is to invert for azimuthal anisotropy and phase velocity simultaneously using all earthquakes at once (**a6_b_eikonal_eq_2DanisoRT.m**). The advantage is that both terms are parameterized explicitly in the inversion. This single-step solution utilizes all data at once so outliers and bad measurements should be less problematic. In addition, for the single-step the smoothing you apply is the smoothing you get in the final map, whereas for the per-event inversion you control only smoothing of each individual event map but the final stacked map can end up much smoother (in my experience). The disadvantage is that the Helmholtz correction cannot be applied. 

There are options for both a 1-D anisotropy (one anisotropy strength and fast azimuth for the entire array) and 2-D anisotropy (anisotropy strength and azimuth at each grid point).

### Attenuation (*WORK IN PROGRESS*)

Rayleigh-wave attenuation is solved for following Bao et al. (2016) GJI [doi: 10.1093/gji/ggw151](https://academic.oup.com/gji/article/206/1/241/2606522) in **b1_estimate_alpha_beta.m**. This requires both traveltime and amplitude fields and their spatial derivatives. We have included the option to apply station corrections following Eddy and Ekström [https://doi.org/10.1016/j.epsl.2014.01.013](https://www.sciencedirect.com/science/article/pii/S0012821X14000223?via%3Dihub) in **a8_a_receiver_terms.m**, under the assumption of an even azimuthal distribution of events.

The calculation requires gradients and Laplacian of the traveltime and amplitude fields. We've included two ways to estimate those fields. The first follows the original ASWMS package technique of fitting a surface to the amplitudes or traveltimes and simply taking the first and second spatial derivatives. However, we have found that this is not always a stable calculation in the presence of noise. Additionally, it is difficult to control the smoothness of the final Laplacian map by smoothing the original surface. 

The second (new) method utilizes the eikonal inversion which effectively solves for the gradient in traveltime. We apply an identical "eikonal" inversion to the amplitude data in order to solve for the gradient in amplitude (**a8_b_ampgrad_inv.m**) and simply take the first spatial derivatives of those maps to get the Laplacian. The advantage of this technique is that the regularized inversion provides a much more stable measurement in the presence of noise. In addition, the smoothness of the Laplacian map can be directly controlled by the 1st derivative smoothing constraint in the inversion.


