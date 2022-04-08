% Run fortran tpw routine

clear;

is_overwrite = 0; % overwrite previous results?

setup_parameters_tpw;
path2bin = parameters_tpw.path2bin;
workingdir_tpw = parameters_tpw.workingdir;
periods = parameters.periods;
currentdir = pwd;

% system(['cp stationid.dat ',workingdir_tpw]);
cd(workingdir_tpw)
for ip = 1:length(periods)
    period = periods(ip);
    
    phvfile = [workingdir_tpw,'/','outvel.',num2str(round(period),'%03d'),'.txt'];
    if exist(phvfile) && ~is_overwrite
        disp(['Already processed ',num2str(period),'s... skipping']);
        continue
    end
    
    % Remove previous files in working directory
    delete([workingdir_tpw,'/*',num2str(round(period),'%03d'),'.inp.sa360kern']);
    
    % Run TPW fortran binary. May need to add permissions via 'chmod ++x simannerr360.gsdf'
    infile = ['all.',num2str(round(period),'%03d')];
%     [stat, log] = system([path2bin,'/srchwave589.JdF.nophase2.iarea < ',infile,' > tpw.log']);
    [stat, ~] = system([path2bin,'/srchwave589.JdF.nophase2.iarea < ',infile]);
    if stat ~= 0     
        error( 'something is wrong at srchwave589... try running a00_make_fortran_executable.m')
    end
end
delete('followit12c');
cd(currentdir)

