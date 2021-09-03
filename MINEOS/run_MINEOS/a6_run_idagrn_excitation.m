% Run idagrn6_sac_excite to calculate full synthetic seismograms and get
% mode excitation for each branch
%
% Must first run run_mineos_check and mk_kernels
%
% JBR 07/18
clear; close all;

COMP = 'T'; %'Z' 'R' 'T' % Component

parameter_FRECHET;
TYPE = param.TYPE;
CARDID = param.CARDID;
EVTPATH = param.EVTPATH;
STAPATH = param.STAPATH;
SYNTH_OUT = param.SYNTH_OUT;

if ( TYPE == 'T') 
    TYPEID = param.TTYPEID;
    if COMP == 'Z' || COMP == 'R'
        error('Toroidal mode: Component must be T');
    end
elseif ( TYPE == 'S') 
    TYPEID = param.STYPEID;
    if COMP == 'T'
        error('Spheroidal mode: Component must be Z or R');
    end
end

%% Change environment variables to deal with gfortran
setenv('GFORTRAN_STDIN_UNIT', '5') 
setenv('GFORTRAN_STDOUT_UNIT', '6') 
setenv('GFORTRAN_STDERR_UNIT', '0')

EXCITE_OUT = [SYNTH_OUT,'excitation'];
if ~exist(EXCITE_OUT)
    mkdir(EXCITE_OUT);
end

%% Run plot_wk
setpath_plotwk;
write_plotwk(TYPE,CARDID)

com = ['cat run_plotwk.',lower(TYPE),' | plot_wk > plot_wk.LOG'];
[status,log] = system(com);
if status ~= 0     
    error( 'something is wrong at plot_wk')
end

%% Run idagrn6_sac_excite
setpath_idagrn;
write_idagrn(TYPE,CARDID,EVTPATH,STAPATH,LENGTH_HR,DT,COMP)

fprintf('------- Calculating Excitation %s0-%d: %s-------\n',TYPE,N_modes-1,COMP)
system(['cat run_idagrn.',lower(TYPE),' > idagrn.in']);
com = ['cat idagrn.in | idagrn6_excite > idagrn.LOG'];
[status,log] = system(com);
if status ~= 0     
    error( 'something is wrong at idagrn6_excite')
end

system(sprintf('mv *.excite.asc %s',EXCITE_OUT));

%% Change the environment variables back to the way they were
setenv('GFORTRAN_STDIN_UNIT', '-1') 
setenv('GFORTRAN_STDOUT_UNIT', '-1') 
setenv('GFORTRAN_STDERR_UNIT', '-1')

delete('idagrn.in','idagrn.LOG','plot_wk.LOG',['*.',lower(TYPE)])