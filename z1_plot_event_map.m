clear;

setup_parameters

per = 55; % period to plot

% input path and files
workingdir = parameters.workingdir;
eventcs_path = [workingdir,'CSmeasure/'];
eikonal_data_path = [workingdir,'eikonal/'];
eikonal_stack_file = [workingdir,'eikonal_stack_',parameters.component];
helmholtz_path = [workingdir,'helmholtz/'];
traveltime_path = [workingdir,'traveltime/'];
periods = parameters.periods;

figdir = [workingdir,'/figs/eventmaps/'];

Iper_plot = find(periods==per);
lalim = parameters.lalim;
lolim = parameters.lolim;

%%

eventfiles = dir([eikonal_data_path,'/*_eikonal_',parameters.component,'.mat']);
evlas=[]; evlos=[]; Mag=[]; baz=[];
ii = 0;
for ie = 1:length(eventfiles)
    load(fullfile(eikonal_data_path,eventfiles(ie).name));
    if isempty(find(~isnan(eventphv(Iper_plot).GV)))
        continue
    end
    ii = ii + 1;
	evlas(ii) = eventphv(Iper_plot).evla;
    evlos(ii) = eventphv(Iper_plot).evlo;
    Mag(ii) = eventphv(Iper_plot).Mw;
    baz(ii) = azimuth(mean(lalim), mean(lolim), eventphv(Iper_plot).evla, eventphv(Iper_plot).evlo);
end

%% Plot map

figure(3); clf;
landareas = shaperead('landareas.shp','UseGeoCoords',true);
ax = axesm('eqdazim', 'Frame', 'on', 'Grid', 'off');
box off;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
% setm(ax,'Origin',[mean(lalim),mean(lolim)])
setm(ax,'Origin',[mean(lalim),mean(lolim)],'FLatLimit',[-125 125]+mean(lalim),'FLonLimit',[],'MapLonLimit',[-125 125]+mean(lolim))
% setm(ax,'Origin',[mean(lalim),mean(lolim)],'FLatLimit',[-125 125]+mean(lalim),'FLonLimit',[-125 125]+mean(lolim))
geoshow(ax, landareas,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); hold on;
plotm(evlas,evlos,'o','color',[0.4 0 0],'MarkerFaceColor',[0.85 0 0],'MarkerSize',10,'linewidth',1);
plotm(mean(lalim),mean(lolim),'p','color',[0 0.2 0.4],'MarkerFaceColor',[0 0.5 1],'MarkerSize',24,'linewidth',1);
for ii = [30 60 90 120]
    [latc,longc] = scircle1(mean(lalim),mean(lolim),ii);
    plotm(latc,longc,'-k','linewidth',1)
end

figure(4); clf;
binc = [0:30:360]*pi/180;
polarhistogram(baz*pi/180,binc,'FaceColor',[0.7 0 0]);
set(gca,'fontsize',16,'linewidth',1.5)
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
thetatickformat('degrees')

disp([num2str(length(evlas)),' events'])

if ~exist(figdir)
    mkdir(figdir)
end
save2pdf([figdir,'events',num2str(periods(Iper_plot)),'s.pdf'],3,500);
save2pdf([figdir,'events_baz_',num2str(periods(Iper_plot)),'s.pdf'],4,500);
