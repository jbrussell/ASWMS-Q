% script to remove events that are close in the time that may interference with each other.
% written by Ge Jin,jinwar@gmail.com


clear;

% Setup parameters
setup_parameters

% eventmatpath = './eventmat/';
% csmatpath = './CSmeasure/';
% eikonalpath = './eikonal/';
workingdir = parameters.workingdir;
eventmatpath = [workingdir,'eventmat/'];
csmatpath = [workingdir,'CSmeasure/'];
eikonalpath = [workingdir,'eikonal/'];

% Setup Error Codes for Bad data
setup_ErrorCode

min_sta_num = parameters.min_sta_num;
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
gridsize = gridsize*3;
xnode = [lalim(1),mean(lalim),lalim(2)];
ynode = [lolim(1),mean(lolim),lolim(2)];
[xi yi] = meshgrid(xnode,ynode);

matfiles = dir([eventmatpath,'/*_',parameters.component,'.mat']);
for ie = 1:length(matfiles)
	% read in the events information
	temp = load([eventmatpath,matfiles(ie).name]);
	event = temp.event;
	disp(event.id)
	evotimes(ie) = event.otime;
	evlas(ie) = event.evla;
	evlos(ie) = event.evlo;
	evids(ie) = {event.id};
	isgood(ie) = 1;
	dist = deg2km(distance(xi,yi,evlas(ie),evlos(ie)));
	win_start(ie) = evotimes(ie) + min(dist(:)/5);
	win_end(ie) = evotimes(ie) + max(dist(:)/2);
	if ~isfield(event,'stadata')
		isgood(ie)=0;
	elseif length(event.stadata)<min_sta_num
		isgood(ie)=0;
	end
end % end of event loop

for ie = 1:length(evlas)
	for je = 1:length(evlas)
		if ie == je
			continue;
		end
		if win_start(ie) > win_start(je) && win_start(ie) < win_end(je)
			isgood(ie) = 0;
			isgood(je) = 0;
		end
		if win_end(ie) > win_start(je) && win_end(ie) < win_end(je)
			isgood(ie) = 0;
			isgood(je) = 0;
		end
	end
end

disp('Bad events:')
badind = find(isgood == 0);
for ie = badind
	disp(evids(ie));
end

%com = input('Do you want to delete these events? y/n','s');
com = 'y';

if com == 'y'
for ie = badind
	delete([eventmatpath,char(evids(ie)),'*.mat'])
	delete([csmatpath ,char(evids(ie)),'*.mat'])
	delete([eikonalpath ,char(evids(ie)),'*.mat'])
end
end

Ngood = length(char(evids(isgood == 1)));
Nevs = length(char(evids));
disp([num2str(Ngood),'/',num2str(Nevs),' events kept']);

%% Plot map

figure(3); clf;
landareas = shaperead('landareas.shp','UseGeoCoords',true);
ax = axesm('eqdazim', 'Frame', 'on', 'Grid', 'off');
box off;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
% setm(ax,'Origin',[mean(lalim),mean(lolim)])
setm(ax,'Origin',[mean(lalim),mean(lolim)],'FLatLimit',[-125 125]+mean(lalim),'FLonLimit',[],'MapLonLimit',[-125 125]+mean(lolim))
geoshow(ax, landareas,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); hold on;
plotm(evlas(logical(isgood)),evlos(logical(isgood)),'o','color',[0.4 0 0],'MarkerFaceColor',[0.85 0 0],'MarkerSize',10,'linewidth',1);
isbad = ~logical(isgood);
if ~isempty(evlas(isbad))
    plotm(evlas(isbad),evlos(isbad),'o','color',[0.4 0 0],'MarkerSize',10,'linewidth',1);
end
plotm(mean(lalim),mean(lolim),'p','color',[0 0.2 0.4],'MarkerFaceColor',[0 0.5 1],'MarkerSize',24,'linewidth',1);
for ii = [30 60 90 120]
    [latc,longc] = scircle1(mean(lalim),mean(lolim),ii);
    plotm(latc,longc,'-k','linewidth',2)
end
