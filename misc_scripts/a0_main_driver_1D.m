%% The driver program to run all the component in this folder
% Written by Ge Jin,jinwar@gmail.com

% if using sac dataset
a1_sac2eventmat

% download the data using 
% data_download

% clean up multiple or close events
a3_cleanup_events

% automatic select the window range
a4_run_autowinpick

% making cross-correlation measurement
a5_gsdfmain

% calculating eikonal tomography for each event
a6_eikonal_eq_flat

% stacking the result
% a7_stack_phv

% apply amplitude correction
% helmholtz_eq

% stack the result of helmholtz
% stack_helm

% export the result into xyz format
% make_xyz

% Calculate and plot 1-D azimuthal anisotropy
average_event_phv_JGR18_2thetafit_eik_bin
