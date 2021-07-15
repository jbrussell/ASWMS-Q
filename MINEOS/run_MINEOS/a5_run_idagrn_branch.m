% Run idagrn6_mask to calculate individual synthetic seismograms for each
% mode branch specified
%
% Must first run run_mineos_check and mk_kernels
%
% JBR 07/18
clear; close all;

branches = [0 1]; % Fundamental -> 0
COMP = 'R'; %'Z' 'R' 'T' % Component

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

for branch = branches
    BRANCH_OUT = [SYNTH_OUT,'br',num2str(branch)];
    if ~exist(BRANCH_OUT)
        mkdir(BRANCH_OUT);
    end
    
    %% Run plot_wk
    setpath_plotwk;
    write_plotwk_branch(TYPE,CARDID,branch)

    com = ['cat run_plotwk_branch.',lower(TYPE),' | plot_wk > plot_wk.LOG'];
    [status,log] = system(com);
    if status ~= 0     
        error( 'something is wrong at plot_wk')
    end

    %% Run idagrn6_mask
    setpath_idagrn;
    write_idagrn_branch(TYPE,CARDID,branch,EVTPATH,STAPATH,LENGTH_HR,DT,COMP)
    
    fprintf('-------Calculating synthetics %s%d: %s-------\n',TYPE,branch,COMP)
    system(['cat run_idagrn_branch.',lower(TYPE),' > idagrn.in']);
    com = ['cat idagrn.in | idagrn6_mask > idagrn.LOG'];
    [status,log] = system(com);
    if status ~= 0     
        error( 'something is wrong at idagrn6_mask')
    end
    

    system(sprintf('mv *.%s.sac %s',COMP,BRANCH_OUT));
end

%% Change the environment variables back to the way they were
setenv('GFORTRAN_STDIN_UNIT', '-1') 
setenv('GFORTRAN_STDOUT_UNIT', '-1') 
setenv('GFORTRAN_STDERR_UNIT', '-1')

delete('idagrn.in','idagrn.LOG','plot_wk.LOG',['*.',lower(TYPE)])