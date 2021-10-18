% Run fortran tpw routine

clear;

setup_parameters_tpw;
path2bin = parameters_tpw.path2bin;
workingdir_tpw = parameters_tpw.workingdir;
periods = parameters.periods;
currentdir = pwd;

% system(['cp stationid.dat ',workingdir_tpw]);
cd(workingdir_tpw)
for ip = 1:length(periods)
    period = periods(ip);
    % Remove previous files in working directory
    delete([workingdir_tpw,'/*',num2str(period,'%03d'),'.inp.sa360kern']);
    
    % Run TPW fortran binary. May need to add permissions via 'chmod ++x simannerr360.gsdf'
    infile = ['all.',num2str(period,'%03d')];
%     [stat, log] = system([path2bin,'/srchwave589.JdF.nophase2.iarea < ',infile,' > tpw.log']);
    [stat, ~] = system([path2bin,'/srchwave589.JdF.nophase2.iarea < ',infile]);
end
cd(currentdir)
