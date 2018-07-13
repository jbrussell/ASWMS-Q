% JBR - 1/7/16
%
% Plot ray geometry before and after quality control step
%
clear;

addpath('../');
setup_parameters;
periods = parameters.periods;
lalim = parameters.lalim;
lolim = parameters.lolim;

ev = 19;

is_QC = 1; % Quality Control?
is_fig = 0;

% --------------------------

eventfile = parameters.eventfile;
evmat = '../eventmat/';
CSmat = '../CSmeasure/';
% fig_PATH = ['./',dir_PROJ,'/figs/'];

% system(['mkdir ',fig_PATH]);

%% LOAD DATA STRUCTURES

% LOAD EVENTS .mat
obs_evpath = [evmat];
obs_evfiles = dir([obs_evpath,'*.mat']);
numev = length(obs_evfiles);

events_obs(numev) = struct();
events_synth(numev) = struct();
for iev = 1:numev
    load([obs_evpath,obs_evfiles(iev).name]); %loads structure called "event"
    events_obs(iev).event = event;
end


% LOAD MEASUREMENTS .mat
obs_CSpath = [CSmat];
obs_CSfiles = dir([obs_CSpath,'*.mat']);
numCS = length(obs_CSfiles);

for iev = 1:numCS
    load([obs_CSpath,obs_CSfiles(iev).name]); %loads structure called "eventcs"
    events_obs(iev).eventcs = eventcs;
end


% QC STEP
if is_QC == 1
events_obs_QC = events_obs;
for iev = 1:numCS
    num_measures = length(events_obs(iev).eventcs.CS);
    for imeas = 1:num_measures
        I_bad_obs = events_obs(iev).eventcs.CS(imeas).isgood ~= 1;
        %I_bad_obs = events_obs(iev).eventcs.CS(imeas).isgood == -2; % high tp error
        events_obs_QC(iev).eventcs.CS(imeas).dtp(I_bad_obs) = NaN;
    end
end
end

%% Plot Ray Geometry

fig20 = figure(20); clf; hold on; set(gcf, 'Color', 'w'); box on;
for imeas = 1:num_measures
    sta1 = events_obs(ev).eventcs.CS(imeas).sta1;
    stla1 = events_obs(ev).event.stadata(sta1).stla;
    stlo1 = events_obs(ev).event.stadata(sta1).stlo;
    sta2 = events_obs(ev).eventcs.CS(imeas).sta2;
    stla2 = events_obs(ev).event.stadata(sta2).stla;
    stlo2 = events_obs(ev).event.stadata(sta2).stlo; 
    plot([stlo1 stlo2],[stla1 stla2],'-k'); hold on;
end
numsta = length(events_obs(ev).event.stadata);
for ista = 1:numsta
    stla = events_obs(ev).event.stadata(ista).stla;
    stlo = events_obs(ev).event.stadata(ista).stlo;
    plot(stlo,stla,'vr','markerfacecolor','r','markeredgecolor','none','markersize',8);
end
title('All rays','fontsize',12);
xlim(lolim);
ylim(lalim);

fig21 = figure(21); clf; set(gcf, 'Color', 'w','position',[145          29        1045         676]);
for iper = 1:length(periods)
    subplot(3,ceil(length(periods)/3),iper); box on; hold on;
    for imeas = 1:num_measures
        if ~isnan(events_obs_QC(ev).eventcs.CS(imeas).dtp(iper))
            sta1 = events_obs(ev).eventcs.CS(imeas).sta1;
            stla1 = events_obs(ev).event.stadata(sta1).stla;
            stlo1 = events_obs(ev).event.stadata(sta1).stlo;
            sta2 = events_obs(ev).eventcs.CS(imeas).sta2;
            stla2 = events_obs(ev).event.stadata(sta2).stla;
            stlo2 = events_obs(ev).event.stadata(sta2).stlo;
            plot([stlo1 stlo2],[stla1 stla2],'-k');
        end
    end
    for ista = 1:numsta
        stla = events_obs(ev).event.stadata(ista).stla;
        stlo = events_obs(ev).event.stadata(ista).stlo;
        plot(stlo,stla,'vr','markerfacecolor','r','markeredgecolor','none','markersize',8);
    end
    xlim(lolim);
    ylim(lalim);
    title([num2str(periods(iper)),' s'],'fontsize',12);
end

if is_fig == 1
    export_fig(fig20,[fig_PATH,'rays_all'],'-pdf','-painters');
    export_fig(fig21,[fig_PATH,'rays_',events_obs(ev).event.id,'_',num2str(periods(per)),'s'],'-pdf','-painters');
end