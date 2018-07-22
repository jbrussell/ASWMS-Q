% JBR - 1/9/16
%
% Plot phase delay (dtp) as a function of epicentral distance difference
% (ddist) similar to Figure 6 in the ASWMS manual.
%

clear;
addpath('../')

setup_parameters;

periods = parameters.periods;

is_fig = 1;
% --------------------------

eventfile = parameters.eventfile;
workingdir = parameters.workingdir;
evmat = ['../',workingdir,'eventmat/'];
CSmat = ['../',workingdir,'CSmeasure/'];
fig_PATH = ['../',workingdir,'figs/dtp_ddist/'];

maxstadist = parameters.maxstadist;
xlims = [-maxstadist maxstadist];

if ~exist(fig_PATH)
    mkdir(fig_PATH);
end
%% LOAD DATA STRUCTURES

% LOAD MEASUREMENTS .mat
obs_CSpath = [CSmat];
obs_CSfiles = dir([obs_CSpath,'*.mat']);
numCS = length(obs_CSfiles);

evmat_files = dir([evmat,'*.mat']);

for iev = 1:numCS
    load([obs_CSpath,obs_CSfiles(iev).name]); %loads structure called "eventcs"
    events_obs(iev).eventcs = eventcs;
    eventmat(iev) = load([evmat,evmat_files(iev).name]);
end

%% Plot phase delays

latarr = mean(parameters.lalim);
lonarr = mean(parameters.lolim);


%Loop over periods for a given event
max_dtp = -1000;
min_dtp = 1000;

% Data
fig56 = figure(56); clf; hold on; set(gcf, 'Color', 'w'); box on;
numper = length(periods);
clr_per = jet(numper);
lgd = {};
h = zeros(numper,1);

for ev = 1:numCS
    azi(ev) = azimuth(latarr,lonarr,eventmat(ev).event.evla,eventmat(ev).event.evlo);    
    avgphv(ev,:) = events_obs(ev).eventcs.avgphv;
end
azi(azi>360) = azi(azi>360) - 360;
azi(azi<0) = azi(azi<0) + 360;

for iper = 1:numper    
    subplot(3,4,iper);
    plot(azi,avgphv(:,iper),'xr');
    xlabel('Azimuth');
    ylabel('phv');
end
% pause;

%% EXPORT FIGURES
if is_fig == 1
%     export_fig(fig55,[fig_PATH,'obs_',num2str(events_obs(ev).eventcs.id)],'-pdf','-painters');
    save2pdf([fig_PATH,'obs_',num2str(events_obs(ev).eventcs.id)],fig56,500)
%     export_fig(fig66,[fig_PATH,'synth_',num2str(events_synth(ev).event.id)],'-pdf','-painters');
end
%pause;