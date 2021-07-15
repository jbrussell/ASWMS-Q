%% Setup Parameters for running MINEOS to calculate senstivity kernels, dispersion, and synthetics

%clear all;
addpath('./functions'); % Path to matlab unctions
path2runMINEOS = './'; % Path to this folder
path2BIN = '../FORTRAN/bin'; % Path to fortran binaries

% Mineos table parameters
maxN = 400000; % Estimate of max number of modes
minF = 0;
maxF = 66.67; % max frequency in mHz; %10.1; %250.05; %333.4; %500.05; %200.05; %%150.05; %50.05;
minL = 0;
maxL = 50000;
N_modes = 10; % <0 uses all mode branches, 1=fundamental only -------- JOSH 8/22/15
param.CARDID = 'pa5_5km'; %'prem_35'; %'fail_H01221_90L'; %'prem_35'; %'Nomelt_taper_aniso_constxicrman_etaPREM_constxilays'; %'pa5_5km';
% param.CARDID = 'prem_35';

% (1 => yes, 0 => no)
SONLY = 1; %Spheroidal modes? (RAYLEIGH)
TONLY = 0; %Toroidal modes? (LOVE)

% for plotting kernels
param.periods = [20, 23, 27, 31, 36, 41, 47, 55, 63, 73, 84];

ch_mode = 0; % (DO NOT CHANGE) mode branch to check for missed eigenfrequencies 0 => T0 ------- JOSH 10/7/15

%% Parameters for idagrn synthetics
LENGTH_HR = 1.0; %1.0; % length of seismogram in hours
DT = 1.0; % 1/samplerate
eventfile = 'evt_201404131236';
stationfile = 'stations.stn';


%%
% Setup idagrn paths
param.IDAGRN = [path2runMINEOS,'/IDAGRN/'];
param.EVTPATH = [param.IDAGRN,'EVT_FILES/',eventfile];
param.STAPATH = [param.IDAGRN,'STATION/',stationfile];
param.SYNTH_OUT = [param.IDAGRN,'SYNTH/',param.CARDID,'_b',num2str(N_modes),'/',eventfile,'/'];
if ~exist(param.SYNTH_OUT)
    mkdir(param.SYNTH_OUT);
end

%%
if SONLY == 1 && TONLY == 0
    param.TYPE = 'S';
elseif SONLY == 0 && TONLY == 1
    param.TYPE = 'T';
else
    error('Choose SONLY or TONLY, not both');
    
end

% Setup Parameters for Initial Model
param.CARD = [param.CARDID,'.card'];
param.CARDPATH  = [path2runMINEOS,'/CARDS/'];
param.TABLEPATH = [path2runMINEOS,'/MODE/TABLES/'];
param.MODEPATH  = [path2runMINEOS,'/MODE/TABLES/MODE.in/'];
if ~exist(param.MODEPATH)
    mkdir(param.MODEPATH);
end
param.RUNPATH = pwd;

%% create dir for output MINEOS automatically, doesn't need to be changed.
CARDTABLE = [param.TABLEPATH,param.CARDID,'/tables/'];
if ~exist(CARDTABLE)
    mkdir([param.TABLEPATH,param.CARDID])
    mkdir(CARDTABLE)
end

%% setup Parameters for kernals
param.frechet = [path2runMINEOS,'/MODE/FRECHET/'];
param.frechetpath = [path2runMINEOS,'/MODE/FRECHET/',param.CARDID,'/'];

if ~exist(param.frechetpath) 
    mkdir(param.frechetpath)
end

%% setup Parameters for eigenfunctions
param.eigpath = [path2runMINEOS,'/MODE/EIGEN/',param.CARDID,'/'];

if ~exist(param.eigpath) 
    mkdir(param.eigpath)
end

%% setup Parameters for Dispersion
param.disperspath = [path2runMINEOS,'/MODE/DISPERSION/',param.CARDID,'/'];

if ~exist(param.disperspath) 
    mkdir(param.disperspath)
end

%% Turn on if only want to calculate S or T or both for mineous
param.SMODEIN = ['s.mode',num2str(floor(minF)),'_',num2str(floor(maxF)),'_b',num2str(N_modes)];
param.STYPEID = ['s',num2str(floor(minF)),'to',num2str(floor(maxF))];
param.TMODEIN = ['t.mode',num2str(floor(minF)),'_',num2str(floor(maxF)),'_b',num2str(N_modes)];
param.TTYPEID = ['t',num2str(floor(minF)),'to',num2str(floor(maxF))];%'t0to150';

%% Setup paths to FORTRAN binaries
PATH = getenv('PATH');
if isempty(strfind(PATH,path2BIN))
%     setenv('PATH', [PATH,':',path2BIN]);
    setenv('PATH', [path2BIN,':',PATH]);
end
