clear

% setup parameters
setup_parameters
setup_ErrorCode

A_Amax_minthresh = 0.6; % Threshold for minimum excitation ratio

workingdir = parameters.workingdir;
receiverterms_path = [workingdir];
% input path
eventcs_path = [workingdir,'CSmeasure/'];
nodal_eventcs_path = [workingdir,'CSmeasure_nodal/'];
comp = parameters.component;
periods = parameters.periods;
lalim = parameters.lalim;
lolim = parameters.lolim;

csmatfiles = dir([eventcs_path,'/*cs_',comp,'.mat']);
ikill = 0;
evlas=[]; evlos=[];
is_node = false(length(csmatfiles),1);
for ie = 1:length(csmatfiles)
	% read in data and set up useful variables
	temp = load([eventcs_path,csmatfiles(ie).name]);
	eventcs =  temp.eventcs;
	disp(eventcs.id)
    
    evlas(ie) = eventcs.evla;
    evlos(ie) = eventcs.evlo;
    
    all_ratios = [];
    ii = 0;
    for ip = 1:length(periods)
        for ista = 1:length(eventcs.autocor)
            ii = ii + 1;
            all_ratios(ip,ista) = eventcs.source(ista).excitation(ip).ratio_AmpMax(1);
        end
    end
    
    for ista = 1:length(eventcs.autocor)
        Ibad_pers = find(all_ratios(:,ista)<A_Amax_minthresh);
        eventcs.autocor(ista).exitflag(Ibad_pers) = ErrorCode.near_node;
        if ~isempty(Ibad_pers)
            ikill = ikill + 1;
            is_node(ie) = true;
            break
        end
    end
    
end % end of loop ie

disp(['Bad: ',num2str(ikill)]);
disp(['Good: ',num2str(ie-ikill)]);

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
plotm(evlas(~is_node),evlos(~is_node),'o','color',[0.4 0 0],'MarkerFaceColor',[0.85 0 0],'MarkerSize',10,'linewidth',1);
if ~isempty(evlas(is_node))
    plotm(evlas(is_node),evlos(is_node),'o','color',[0.4 0 0],'MarkerSize',10,'linewidth',1);
end
plotm(mean(lalim),mean(lolim),'p','color',[0 0.2 0.4],'MarkerFaceColor',[0 0.5 1],'MarkerSize',24,'linewidth',1);
for ii = [30 60 90 120]
    [latc,longc] = scircle1(mean(lalim),mean(lolim),ii);
    plotm(latc,longc,'-k','linewidth',2)
end