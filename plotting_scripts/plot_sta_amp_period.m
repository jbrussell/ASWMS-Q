% code to plot the average amplitude of all stations

clear
addpath('../');addpath('../functions/');
setup_parameters;
comp = parameters.component;
periods = parameters.periods;

avg_band = 1:11; %12:15; %3:6;

outfile = 'NoMelt_sta_info';
outfile_corr = 'NoMelt_sta_info_corr';
station_file = 'stalist.txt';

% CSfiles = dir(['CSmeasure/*_',comp,'*.mat']);
workingdir = ['../',parameters.workingdir];
CSfiles = dir([workingdir,'CSmeasure/*_',comp,'*.mat']);

% Gather information
stnms = {};
stainfo = [];
for ie = 1:length(CSfiles)
	clear eventcs amps
	amps = [];
% 	load(fullfile('CSmeasure',CSfiles(ie).name));
    load(fullfile(workingdir,'CSmeasure',CSfiles(ie).name));
	disp(CSfiles(ie).name)
    icount = 0;
	for ista = 1:length(eventcs.stnms)
        if isempty(find(eventcs.autocor(ista).exitflag>0))
            continue
        end
        icount = icount+1;
        amps(icount,:) = sqrt(eventcs.autocor(ista).amp);
        % Remove geometrical spreading effect on amplitude
        geomsprd = sqrt(abs(sind(km2deg(eventcs.dists(ista)))));
        amps(icount,:) = amps(icount,:) * geomsprd;
	end
	meanamp = mean(amps,1);
    
    icount = 0;
	for ista = 1:length(eventcs.stnms)
        if isempty(find(eventcs.autocor(ista).exitflag>0))
            continue
        end
        icount = icount+1;
		if ismember(eventcs.stnms(ista),stnms)
			staid = find(ismember(stnms,eventcs.stnms(ista)));
			stainfo(staid).norm_amp(end+1,:) = amps(icount,:)./meanamp;
			stainfo(staid).ori_amp(end+1,:) = amps(icount,:);
            stainfo(staid).stlas = eventcs.stlas(ista);
            stainfo(staid).stlos = eventcs.stlos(ista);
		else
			stnms(end+1) = eventcs.stnms(ista);
			stainfo(end+1).stnm = eventcs.stnms(ista);
			stainfo(end).norm_amp(1,:) = amps(icount,:)./meanamp;
			stainfo(end).ori_amp(1,:) = amps(icount,:);
            stainfo(end).stlas = eventcs.stlas(ista);
            stainfo(end).stlos = eventcs.stlos(ista);
		end
	end
end

save([outfile,'.mat'],'stainfo');

load(outfile)

% badstnms = textread('badampsta.lst','%s');
% OBSstnms = textread(station_file,'%s');

% calculate the means
for ista = 1:length(stainfo)
	stainfo(ista).meanamp = mean(stainfo(ista).norm_amp,1);
	stainfo(ista).avgmean = mean(stainfo(ista).meanamp(avg_band));
    if exist('badstnms','var') && ismember(stainfo(ista).stnm,badstnms)
        stainfo(ista).isgood = 0;
    else
        stainfo(ista).isgood = 1;
    end
    if exist('OBSstnms','var') && ismember(stainfo(ista).stnm,OBSstnms)
        stainfo(ista).isOBS = 1;
    else
        stainfo(ista).isOBS = 0;
    end
end

%calculate the correction
isgood = [stainfo.isgood];
isOBS = [stainfo.isOBS];
onland_ind = find(isgood==1 & isOBS == 0);
onland_avg = mean([stainfo(onland_ind).avgmean]);
OBS_ind = find(isgood==1 & isOBS ==1);
OBS_avg = mean([stainfo(OBS_ind).avgmean]);
OBS_cor = onland_avg/OBS_avg;

%% making plots
% Before the correction
figure(23)
clf
subplot(2,1,1)
hold on
N = ceil(sqrt(length(periods)));
M = ceil(sqrt(length(periods)));
for ip = 1:length(periods)
    subplot(N,M,ip); hold on;
	for ista = 1:length(stainfo)
		norm_amp = stainfo(ista).norm_amp;
		x = ones(size(norm_amp,1),1)*ista;
		plot(x,norm_amp(:,ip),'x','color',[0.8 0.8 0.8]);
		errorbar(ista, mean(norm_amp(:,ip)), std(norm_amp(:,ip)),'ro','markerfacecolor','r');
		if stainfo(ista).isgood
			plot(ista, stainfo(ista).avgmean,'bo','markerfacecolor','b','markersize',3);
		end
		maxy = 3;
		ylim([0 maxy]);
%		text(ista,-maxy/10,char(stainfo(ista).stnm),'rotation',90);
		text(ista,-maxy/50,char(stainfo(ista).stnm),'rotation',-90);
		set(gca, 'XTickLabel','')
        title([num2str(periods(ip)),' s']);
	end
end
% subplot(N,M,ip+1); hold on;
% for ista = 1:length(stainfo)
%     if stainfo(ista).isgood
%         plot(ista, stainfo(ista).avgmean,'go','markerfacecolor','g','markersize',5);
%     end
%     ylim([0 maxy]);
%     text(ista,-maxy/50,char(stainfo(ista).stnm),'rotation',-90);
%     set(gca, 'XTickLabel','')
%     title('Average');
% end
% subplot(2,1,2)
% isgood = [stainfo.isgood];
% ind = find(isgood);
% hist([stainfo(ind).avgmean],10);
% xlim([0 2])
