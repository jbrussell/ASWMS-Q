% Script to print out grid nodes needed for fortran two plane wave method
% NJA, LDEO 2013
%
% Format is as follows and moves from left to right
% ncol x nrow
% # of grid nodes
% latitude longitude
% corner points x4
% ncol
% col_spacing row_spacing (in km)
setup_parameters_tpw;

% Parameters
workingdir_tpw = parameters_tpw.workingdir;
path2bin = parameters_tpw.path2bin;
refvs = parameters_tpw.refphv;
periods = parameters.periods;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters_tpw.gridsize;
gridsize_interp = parameters_tpw.gridsize_interp;
gridla = lalim(1):gridsize:lalim(2);
gridlo = lolim(1):gridsize:lolim(2);

if gridsize_interp > gridsize
    error('Interpolated grid spacing is larger than inversion grid...');
end

currentdir = pwd;
cd(workingdir_tpw)

%% Build input file for building inversion grid within main TPW code
% This grid is assumed to be that same as that defined in ASWMS
%

slat = mean(gridla); % center lat
slon = mean(gridlo); % center lon
sazimz = 0; % azimuth from center point that will tilt grid relative to North
nxpt = length(gridla); % number of points in latitude
delx = mean(abs(diff(gridla))); % spacing in latitude
begx = min(gridla) - slat; % beginning latitude distance from center
nypt = length(gridlo); % number of points in longitude
dely = mean(abs(diff(gridlo))); % spacing in longitude
begy = min(gridlo) - slon; % beginning longitude distance from center
toplecrn = [max(gridla) min(gridlo)]; % top left corner of region
topricrn = [max(gridla) max(gridlo)]; % top right corner of region
botlecrn = [min(gridla) min(gridlo)]; % bottom left corner of region
botricrn = [min(gridla) max(gridlo)]; % bottom right corner of region

gridinpfile = ['gridinp.dat'];
if exist(gridinpfile)
    delete(gridinpfile);
end
fid = fopen(gridinpfile,'w');
fprintf(fid,'%.2f %.2f\n',slat,slon);
fprintf(fid,'%.2f\n',sazimz);
fprintf(fid,'%d %.2f %.2f\n',nxpt,delx,begx);
fprintf(fid,'%d %.2f %.2f\n',nypt,dely,begy);
fprintf(fid,'%.2f %.2f\n',toplecrn); % top left corner
fprintf(fid,'%.2f %.2f\n',topricrn); % top right corner
% fprintf(fid,'%.2f %.2f\n',botlecrn); % bottom left corner
% fprintf(fid,'%.2f %.2f\n',botricrn); % bottom right corner
fprintf(fid,'%.2f %.2f\n',botricrn); % bottom right corner
fprintf(fid,'%.2f %.2f\n',botlecrn); % bottom left corner

fclose(fid);

%% Generate output node locations and starting model
% This is an interpolated version of the inversion grid. This can be finer
% scale than the inversion grid. Must extend slightly beyond the inversion grid.
dlat = gridsize_interp;
beglat = min(lalim)-dlat*4;
endlat = max(lalim)+dlat*4;
dlon = gridsize_interp;
beglon = min(lolim)-dlon*4;
endlon = max(lolim)+dlon*4;

for ip = 1:length(periods)
    period = periods(ip);
    fname = ['outgrd.',num2str(round(period),'%03d'),'.inp'];
    
    % Generate input file for poutgrd
    infile_poutgrd = [workingdir_tpw,'/poutgrd.in'];
    if exist(infile_poutgrd)
        delete(infile_poutgrd);
    end
    fid = fopen(infile_poutgrd,'w');
    if (fid == -1)
        error (['    Cannot open file: ', infile_poutgrd]);
    end
    fprintf(fid,'%.2f %.2f %.2f %.2f %.2f %.2f\n',beglat,endlat,dlat,beglon,endlon,dlon);
    fprintf(fid,'%.3f\n',refvs(ip));
    fprintf(fid,'%s',fname);
    fclose(fid);

    % Run poutgrd
    [stat, log] = system([path2bin,'/poutgrd < ',infile_poutgrd]);
    if stat ~= 0     
        error( 'something is wrong at poutgrd... try running a00_make_fortran_executable.m')
    end
end

%%
cd(currentdir)

