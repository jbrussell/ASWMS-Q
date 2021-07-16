% OPTIONAL
%
% Script to calculate Rayleigh-wave source excitation using equation
% (11.34) from Dahlen & Tromp. This requires eigenfunctions and dispersion
% for a desired earth model, which must be precalculated using MINEOS. This
% also requires the moment tensor components, which will have been
% saved if you used https://github.com/jbrussell/fetch_EVENTS to
% retrieve the waveform data.
%
% Credit to Anant Hariharan for coding up the excitation equations
%
% github.com/jbrussell
% 7/2021


clear;
addpath('./functions_excitation/');

CARD = 'pa5_5km'; % name of mineos card file
path2eig = ['./MINEOS/run_MINEOS/MODE/EIGEN/',CARD,'/']; % path to eigenfunction location
path2phv = ['./MINEOS/run_MINEOS/MODE/TABLES/',CARD,'/tables/']; % path to mode table location
% Path to idagrn CMT file
path2cmt = '~/BROWN/RESEARCH/PROJ_NoMelt/DATA/EVENTS/fetch_EVENTS/IRIS_ZA_5.5_Zcorr/CMT2idagrn/';
MaxN = 1; % maximum overtone for estimating overtone interference
MinorOrMajor = 0; % Set the variable to 1 for major, 0 for minor arc overtone interference

% Setup parameters
setup_parameters

workingdir = parameters.workingdir;
eventmatpath = [workingdir,'eventmat/'];
csmatpath = [workingdir,'CSmeasure/'];
eikonalpath = [workingdir,'eikonal/'];
periods = parameters.periods;

% Setup Error Codes for Bad data
setup_ErrorCode

lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
gridsize = gridsize*3;
xnode = [lalim(1),mean(lalim),lalim(2)];
ynode = [lolim(1),mean(lolim),lolim(2)];
[xi yi] = meshgrid(xnode,ynode);

matfiles = dir([eventmatpath,'/*_',parameters.component,'.mat']);
for ie = 1:length(matfiles)
	% read in the events information
	temp = load([eventmatpath,matfiles(ie).name]);
	event = temp.event;
    evid = event.id;
	disp(evid)
    
    % Load moment tensor values from CMT idagrn file
    fname = [path2cmt,'/evt_',evid];    
    fid = fopen(fname,'r');
    evid = sscanf(fgetl(fid),'%s');
    temp = sscanf(fgetl(fid),'%f %f %f');
    lat = temp(1);
    lon = temp(2);
    depth_km = temp(3);
    mult_fac = sscanf(fgetl(fid),'%f');
    m_rr = sscanf(fgetl(fid),'%f');
    m_tt = sscanf(fgetl(fid),'%f');
    m_pp = sscanf(fgetl(fid),'%f');
    m_rt = sscanf(fgetl(fid),'%f');
    m_rp = sscanf(fgetl(fid),'%f');
    m_tp = sscanf(fgetl(fid),'%f');    
    fclose(fid);
    
    event.moment.m_rr = m_rr;
    event.moment.m_tt = m_tt;
    event.moment.m_pp = m_pp;
    event.moment.m_rt = m_rt;
    event.moment.m_rp = m_rp;
    event.moment.m_tp = m_tp;
    event.moment.mult_fac = mult_fac;
    
    % Calculate excitation ratios and attach to event structure
    [event] = Rayleigh_GetExcitationRatios(event,periods,CARD,path2phv,path2eig,MaxN,MinorOrMajor);
    
    save([eventmatpath,matfiles(ie).name],'event');
    
end % end of event loop

