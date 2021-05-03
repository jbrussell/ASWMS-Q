% Apply amplitude correction on the result of eikonal tomography
% also use the stacking result for more accurate correction
% written by Ge Jin, jinwar@gmail.com
% 2013-03

clear;

isfigure = 1;
isoverwrite = 1;
is_save_amp_fig = 1;

% setup parameters
setup_parameters

r = 0.05;

% input path and files
% eventcs_path = './CSmeasure/';
% eikonal_data_path = './eikonal/';
% eikonal_stack_file = ['eikonal_stack_',parameters.component];
% attenuation_path = './attenuation/';
workingdir = parameters.workingdir;
eventcs_path = [workingdir,'CSmeasure/'];
eikonal_data_path = [workingdir,'eikonal/'];
eikonal_stack_file = [workingdir,'eikonal_stack_',parameters.component];
helmholtz_path = [workingdir,'helmholtz/'];
attenuation_path = [workingdir,'attenuation/'];

if ~exist(attenuation_path,'dir')
	mkdir(attenuation_path);
end

% load stacked phase velocity map
load(eikonal_stack_file);

% set up useful variables
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
xnode = lalim(1):gridsize:lalim(2);
ynode = lolim(1):gridsize:lolim(2);
[xi yi] = ndgrid(xnode,ynode);
tp_var_tol = parameters.tp_var_tol;
alpha_range = parameters.alpha_range;
alpha_search_grid = parameters.alpha_search_grid;
periods = parameters.periods;

eventfiles = dir([eikonal_data_path,'/*_eikonal_',parameters.component,'.mat']);

load seiscmap

if exist('badampsta.lst','file')
	badstnms = textread('badampsta.lst','%s');
	disp('Found Bad amplitude stations:')
	for ista = 1:length(badstnms)
	disp(badstnms(ista))
	end
end

for ie = 1:length(eventfiles)
%for ie = 59
	% read in data for this event
	clear eventphv eventcs attenuation;
	load(fullfile(eikonal_data_path,eventfiles(ie).name));
	eventid = eventphv(1).id;
	matfilename = fullfile(attenuation_path,[eventphv(1).id,'_attenuation_',parameters.component,'.mat']);
	if exist(matfilename,'file') && ~isoverwrite
		disp(['exist: ',matfilename,', skip!'])
		continue;
	end
	disp(eventid);
	eventcsfile = [eventcs_path,'/',eventid,'_cs_',parameters.component,'.mat'];
	if exist(eventcsfile,'file')
		load(eventcsfile);
	else
		disp(['Cannot find CS file for ',eventid,', Skipped']);
		continue;
    end
    eventhelmholtzfile = [helmholtz_path,'/',eventid,'_helmholtz_',parameters.component,'.mat'];
	if exist(eventhelmholtzfile,'file')
		load(eventhelmholtzfile);
	else
		disp(['Cannot find Helmholtz file for ',eventid,', Skipped']);
		continue;
	end
	if length(eventphv) ~= length(eventcs.avgphv)
		disp('Inconsist of period number for CS file and eikonal file');
		continue;
	end

	for ip = 1:length(eventphv)
		%% fit the amplitude surface
		% reset the arrays
		clear stlas stlos tp
		stlas = eventcs.stlas;
		stlos = eventcs.stlos;
		stnms = eventcs.stnms;
		if exist('badstnms','var')
			list_badstaids = find(ismember(eventcs.stnms,badstnms));
		else
			list_badstaids = [];
		end
		tp = zeros(1,length(stlas));
        if length(eventcs.autocor) ~= length(eventphv(ip).traveltime)
            error('something is wrong... check number of stations');
        end
		for ista = 1:length(eventcs.autocor)
			if eventcs.autocor(ista).exitflag(ip)>0
                tp(ista) = eventphv(ip).traveltime(ista);
			else
				tp(ista) = NaN;
			end
        end

		% get rid of bad stations
		badstaids = find(isnan(tp));
		stlas(badstaids) = [];
		stlos(badstaids) = [];
		tp(badstaids) = [];
		badstanum = 0; badstaids = [];
		for ista = 1:length(tp)
			if stlas(ista) < lalim(1) || stlas(ista) > lalim(2) || ...
					stlos(ista) < lolim(1) || stlos(ista) > lolim(2) || ismember(ista,list_badstaids);
				badstanum = badstanum+1;
				badstaids(badstanum) = ista;
				continue;
			end
			dist = distance(stlas(ista),stlos(ista),stlas,stlos);
			dist = deg2km(dist);
			nearstaids = find(dist > parameters.minstadist & dist < parameters.maxstadist );
			nearstaids(find(ismember(nearstaids,badstaids))) = [];
			if isempty(nearstaids)
				badstanum = badstanum+1;
				badstaids(badstanum) = ista;
				continue;
			end
			meantp = median(tp(nearstaids));
% 			if tp(ista) < meantp./tp_var_tol | tp(ista) > meantp.*tp_var_tol
% 				badstanum = badstanum+1;
% 				badstaids(badstanum) = ista;
% 			end
		end
		stlas(badstaids) = [];
		stlos(badstaids) = [];
		tp(badstaids) = [];
		[tpmap,mesh_xi,mesh_yi]=gridfit_jg(stlas,stlos,tp,xnode,ynode,...
							'smooth',2,'regularizer','del4','solver','normal');

		%% Calculate the traveltime and amplitude fields
		tp_lap=del2m(mesh_xi,mesh_yi,tpmap);
        tp_grad=delm(mesh_xi,mesh_yi,tpmap);
        
        ampmap = helmholtz(ip).ampmap;
        amp_grad=delm(mesh_xi,mesh_yi,ampmap);
        
		% prepare the avg phase velocity and event phase velocity
		avgGV = avgphv(ip).GV;
		if sum(size(avgGV)==size(xi)) < 2
			avgGV = interp2(avgphv(ip).xi,avgphv(ip).yi,avgphv(ip).GV,xi,yi,'linear',NaN);
		end
		eventGV = eventphv(ip).GV;

		% fill in informations
		attenuation(ip).evla = eventphv(ip).evla;
		attenuation(ip).evlo = eventphv(ip).evlo;
        attenuation(ip).Mw = eventphv(ip).Mw;
		attenuation(ip).raydense = eventphv(ip).raydense;
		attenuation(ip).goodnum = eventphv(ip).goodnum;
		attenuation(ip).badnum = eventphv(ip).badnum;
		attenuation(ip).id = eventphv(ip).id;
		attenuation(ip).xi = xi;
		attenuation(ip).yi = yi;
		attenuation(ip).GV_cor = helmholtz(ip).GV_cor;
		attenuation(ip).GV = eventGV;
% 		attenuation(ip).alpha_errs = alpha_errs;
% 		attenuation(ip).alphas = alphas;
% 		attenuation(ip).bestalpha = bestalpha;
        attenuation(ip).ampmap = helmholtz(ip).ampmap;
        attenuation(ip).amp_grad = amp_grad;
        attenuation(ip).amps = helmholtz(ip).amps;
		attenuation(ip).tpmap = tpmap;
        attenuation(ip).tp_lap = tp_lap;
        attenuation(ip).tp_grad = tp_grad;
        attenuation(ip).tp = tp;
		attenuation(ip).period = periods(ip);
% 		bestalphas(ip,ie) = bestalpha;

		% plot to check
		if isfigure
			figure(37)
			clf
                        set(gcf,'renderer','zbuffer');
% 			subplot(2,2,1)
% 			ax = worldmap(lalim, lolim);
% 			surfacem(xi,yi,eventGV);
%             if ~isnan(nanmean(eventGV(:)))
%                 caxis([nanmean(eventGV(:))*(1-r) nanmean(eventGV(:))*(1+r)])
%             end
% 			colorbar
% 			title('before cor');
% 			subplot(2,2,2)
% 			ax = worldmap(lalim, lolim);
% 			surfacem(xi,yi,GV_cor);
%             if ~isnan(nanmean(GV_cor(:)))
%                 caxis([nanmean(GV_cor(:))*(1-r) nanmean(GV_cor(:))*(1+r)])
%             end
% 			colorbar
% 			title('after cor');
			nanind = find(isnan(eventGV(:)));
			tpmap = tpmap';
			tpmap(nanind) = NaN;
            tp_grad = tp_grad';
			tp_grad(nanind) = NaN;
			tp_lap = tp_lap';
			tp_lap(nanind) = NaN;
			subplot(2,2,1)
			ax = worldmap(lalim, lolim);
			surfacem(xi,yi,tpmap);
			title('travel time map')
			plotm(stlas,stlos,'v')
            la_gc = [];
            lo_gc = [];
            for ista = 1:length(stlas)
                [la,lo]=track2('gc',eventphv(ip).evla,eventphv(ip).evlo,stlas(ista),stlos(ista));
                la_gc = [la_gc; la; nan];
                lo_gc = [lo_gc; lo; nan];
            end
            plotm(la_gc,lo_gc,'-k');
            colormap(seiscmap)
			colorbar
            subplot(2,2,3)
			ax = worldmap(lalim, lolim);
			surfacem(xi,yi,tp_grad);
			colorbar
% 			[temp bestalphai] = min(alpha_errs);
			title('\nabla \tau_p')
            drawnow;
			subplot(2,2,4)
			ax = worldmap(lalim, lolim);
			surfacem(xi,yi,tp_lap);
			colorbar
% 			[temp bestalphai] = min(alpha_errs);
			title('\nabla^2 \tau_p')
            drawnow;
           
            
            if is_save_amp_fig
                figure(39);
                if ip == 1
                    clf;
                    set(gcf,'Position',[84           3         744        1022]);
                    sgtitle([eventphv(ip).id,' M',num2str(eventphv(ip).Mw)],'fontweight','bold','fontsize',18)
                    
                    axes('Position',[.4 .005 .35*.6 .4*.6])
                    landareas = shaperead('landareas.shp','UseGeoCoords',true);
                    ax = axesm('eqdazim', 'Frame', 'on', 'Grid', 'off');
                    ax.XAxis.Visible = 'off';
                    ax.YAxis.Visible = 'off';
                    box off;
                    % setm(ax,'Origin',[mean(lalim),mean(lolim)])
                    setm(ax,'Origin',[mean(lalim),mean(lolim)],'FLatLimit',[-125 125]+mean(lalim),'FLonLimit',[],'MapLonLimit',[-125 125]+mean(lolim))
                    geoshow(ax, landareas,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); hold on;
                    for ii = [30 60 90 120]
                        [latc,longc] = scircle1(mean(lalim),mean(lolim),ii);
                        plotm(latc,longc,'-','color',[0.6 0.6 0.6],'linewidth',1)
                    end
                    [la_gcev,lo_gcev]=track2('gc',eventphv(ip).evla,eventphv(ip).evlo,mean(lalim),mean(lolim));
                    plotm(la_gcev,lo_gcev,'-k','linewidth',2);
                    plotm(mean(lalim),mean(lolim),'p','color',[0 0.2 0.4],'MarkerFaceColor',[0 0.5 1],'MarkerSize',24,'linewidth',1);
                    plotm(eventphv(ip).evla,eventphv(ip).evlo,'o','color',[0.4 0 0],'MarkerFaceColor',[0.85 0 0],'MarkerSize',10,'linewidth',1);
                end
                N=3; M = floor(length(periods)/N)+1;
                subplot(M,N,ip)
                ax = worldmap(lalim, lolim);
                surfacem(xi,yi,tpmap);
                plotm(stlas,stlos,'v');
                plotm(la_gc,lo_gc,'-k');
                title([num2str(periods(ip)),' s'],'fontsize',15)
                cb = colorbar;
                clim = cb.Limits;
                colormap(seiscmap)
                
                figure(40);
                if ip == 1
                    clf;
                    set(gcf,'Position',[84           3         744        1022]);
                    sgtitle([eventphv(ip).id,' M',num2str(eventphv(ip).Mw)],'fontweight','bold','fontsize',18)
                    
                    axes('Position',[.4 .005 .35*.6 .4*.6])
                    landareas = shaperead('landareas.shp','UseGeoCoords',true);
                    ax = axesm('eqdazim', 'Frame', 'on', 'Grid', 'off');
                    ax.XAxis.Visible = 'off';
                    ax.YAxis.Visible = 'off';
                    box off;
                    % setm(ax,'Origin',[mean(lalim),mean(lolim)])
                    setm(ax,'Origin',[mean(lalim),mean(lolim)],'FLatLimit',[-125 125]+mean(lalim),'FLonLimit',[],'MapLonLimit',[-125 125]+mean(lolim))
                    geoshow(ax, landareas,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); hold on;
                    for ii = [30 60 90 120]
                        [latc,longc] = scircle1(mean(lalim),mean(lolim),ii);
                        plotm(latc,longc,'-','color',[0.6 0.6 0.6],'linewidth',1)
                    end
                    [la_gcev,lo_gcev]=track2('gc',eventphv(ip).evla,eventphv(ip).evlo,mean(lalim),mean(lolim));
                    plotm(la_gcev,lo_gcev,'-k','linewidth',2);
                    plotm(mean(lalim),mean(lolim),'p','color',[0 0.2 0.4],'MarkerFaceColor',[0 0.5 1],'MarkerSize',24,'linewidth',1);
                    plotm(eventphv(ip).evla,eventphv(ip).evlo,'o','color',[0.4 0 0],'MarkerFaceColor',[0.85 0 0],'MarkerSize',10,'linewidth',1);
                end
                subplot(M,N,ip)
                ax = worldmap(lalim, lolim);
                scatterm(stlas,stlos,100,tp,'v','filled','markeredgecolor',[0 0 0]);
                plotm(la_gc,lo_gc,'-k');
                title([num2str(periods(ip)),' s'],'fontsize',15)
                colorbar
                caxis(clim);
                colormap(seiscmap)
            end
            
%             pause;
		end % end of isfigure
	end  % loop of period
	matfilename = fullfile(attenuation_path,[eventphv(1).id,'_attenuation_',parameters.component,'.mat']);
	save(matfilename,'attenuation');
	fprintf('\n');
	disp(['Saved to ',matfilename]);
    if is_save_amp_fig
        figdir = [workingdir,'/figs/attenuation/'];
        if ~exist(figdir)
            mkdir(figdir);
        end
        save2pdf([figdir,eventphv(1).id,'_attenuation_',parameters.component,'.pdf'],39,100);
        save2pdf([figdir,eventphv(1).id,'_attenuation_',parameters.component,'_StaAmps.pdf'],40,100);
    end
end  % loop of events
