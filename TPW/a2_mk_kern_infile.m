% Script to run seizmo script makekernels to make amplitude and phase
% sensitivity kernels for Yang and Forsyth codes.
% NJA, Summer 2013
clear all

setup_parameters_tpw;

% Parameters
workingdir_tpw = parameters_tpw.workingdir;
width = parameters_tpw.width;
offset = parameters_tpw.offset;
kern_grid_km = parameters_tpw.kern_grid_km;
periods = parameters.periods;
refvs = parameters_tpw.refphv; %[3.53 3.61 3.7 3.76 3.82 3.85];
lalim = parameters.lalim;
lolim = parameters.lolim;

for ip = 1:length(periods)
	period = periods(ip);

	range = [1/period 1/period];

	% First, get some frequency bands:
	band=filter_bank(range,'variable',width,offset);

	% Third, use the beat length as a window width:
	L_beat=1./(band(:,3)-band(:,2));

	% Fourth, modify that window width by an empirically derived
	% power law to get a more typical window width:
	win=L_beat'*(2.5+1000*band(:,1).^2);
    swin = win*[0 1];
    
    % Define time window
    twin = [-2000 7000];
    
    % Get zpad
    zpad=[swin(:,1)-twin(:,1)  twin(:,2)-swin(:,2)];
    
    % setup the grid
    spacing_km = kern_grid_km;
%     width_km = deg2km(max([abs(diff(lalim)) abs(diff(lolim))]))+spacing_km*3; %1000;
%     width_km = floor(width_km/spacing_km)*spacing_km;
    width_km = 6000;
    x=spacing_km:spacing_km:width_km/2;
    x=[-x(end:-1:1) 0 x];
    
    % Get frequency-amplitude values for a windowed Rayleigh wave:
%     [f,a]=getmainlobe(f0,fs,swin,tprfrac,zpad)
    [f,a]=getmainlobe(band(:,1),1,swin,0.2,zpad);
    % Example for 100 s wave
%     [f,a]=getmainlobe(1/100,1,[1000 2000],200/1000,[1000 1000]);
    Nf = length(f);

    % Write out kernel file
    nxkern = length(x);
    xbegkern = x(1);
    dxkern = spacing_km;
    kernelname = [workingdir_tpw,'/sensspec',num2str(period,'%03d'),'s_',num2str(max(x)),'km.dat'];
    
    fid = fopen(kernelname,'w');
    fprintf(fid,'%d %.2f %.2f\n',nxkern,xbegkern,dxkern);
    fprintf(fid,'%d %.2f %.2f\n',nxkern,xbegkern,dxkern);
    fprintf(fid,'%d\n',Nf);
    for ii = 1:Nf
        fprintf(fid,'%d %f\n',f(ii),a(ii));
    end
    fclose(fid);
    
 end
