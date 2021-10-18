% Setup parameters file for two-plane wave

path2ASWMS_output = './ASWMS_OUT/'; % location of ASWMS output directory
path2ASWMS_functions = '../functions'; % location of ASWMS functions
addpath(path2ASWMS_output); addpath(path2ASWMS_functions); addpath('./functions/');
setup_parameters;
parameters.workingdir = path2ASWMS_output;

% Path to tpw fortran binary
parameters_tpw.path2bin = [pwd,'/fortran/'];

% Project directory
parameters_tpw.workingdir = [pwd,'/OUT/'];
if ~exist(parameters_tpw.workingdir)
    mkdir(parameters_tpw.workingdir);
end

% Get reference velocities from Helmholtz stack
temp = load([path2ASWMS_output,'/helmholtz_stack_',parameters.component,'.mat']);
for ip = 1:length(parameters.periods)
    parameters_tpw.refphv(ip) = nanmean(temp.avgphv(ip).GV_cor(:)); % used as the homogeneous starting model
end

% Station file
parameters_tpw.stalist = 'stanumlist';

% Parameters for grid spacing
parameters_tpw.gridsize = parameters.gridsize; % [deg] node spacing for inversion grid
parameters_tpw.gridsize_interp = 0.25; % [deg] node spacing for interpolated output grid

% Parameters for building kernel input file: a2_mk_kern
parameters_tpw.kern_grid_km = 20; % kernel grid spacing in km
% parameters_tpw.smlength = 40; % currently unused as kernels are now calculated within the main TPW fortran routine
parameters_tpw.width = 0.10;
parameters_tpw.offset = 0.1;

% Parameters for conversion from ASWMS to TPW: a3_eikonal_to_TPW
parameters_tpw.min_sta_num = parameters.min_sta_num; % minimum number of station to consider event
parameters_tpw.nfreq = 1; % (DO NOT CHANGE) number of frequencies to compute TPW at a time
parameters_tpw.iterlimit = 10; % maximum # of iterations
parameters_tpw.dampvel = 0.25; % a priori stddev for velocity terms
parameters_tpw.dampaniso = 0.25; % a priori stddev for aniso coefficients
parameters_tpw.divfac = 1.0; % smaller value increases velocity smoothing coefficient
parameters_tpw.divfac_azi = 1.0; % smaller value increases anisotropy smoothing coefficient
