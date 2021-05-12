% Apply amplitude correction on the result of eikonal tomography
% also use the stacking result for more accurate correction
% written by Ge Jin, jinwar@gmail.com
% 2013-03

clear;

isfigure = 1;
isoverwrite = 1;
is_save_amp_fig = 1;

is_receiver_terms = 0; % Use receiver terms calculated from a8_0_receiver_terms ?
is_eikonal_ampgrad = 1; % 1: use eikonal tomography values for amplitude gradient; 0: use amplitude field estimates

% setup parameters
setup_parameters

r = 0.05;

% input path and files
workingdir = parameters.workingdir;
eventcs_path = [workingdir,'CSmeasure/'];
eikonal_data_path = [workingdir,'eikonal/'];
eikonal_stack_file = [workingdir,'eikonal_stack_',parameters.component];
helmholtz_path = [workingdir,'helmholtz/'];
receiverterms_path = [workingdir];
ampgrad_output_path = [workingdir,'ampgrad/'];

if ~exist(helmholtz_path,'dir')
	mkdir(helmholtz_path);
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
amp_var_tol = parameters.amp_var_tol;
alpha_range = parameters.alpha_range;
alpha_search_grid = parameters.alpha_search_grid;
periods = parameters.periods;
min_sta_num = parameters.min_sta_num;


eventfiles = dir([eikonal_data_path,'/*_eikonal_',parameters.component,'.mat']);

load seiscmap

if exist('badampsta.lst','file')
	badstnms = textread('badampsta.lst','%s');
	disp('Found Bad amplitude stations:')
	for ista = 1:length(badstnms)
	disp(badstnms(ista))
	end
end

if is_receiver_terms==1
    load([receiverterms_path,'receiver_terms_',parameters.component,'.mat']);
end

for ie = 1:length(eventfiles)
%for ie = 59
	% read in data for this event
	clear eventphv eventcs helmholtz;
	load(fullfile(eikonal_data_path,eventfiles(ie).name));
	eventid = eventphv(1).id;
	matfilename = fullfile(helmholtz_path,[eventphv(1).id,'_helmholtz_',parameters.component,'.mat']);
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
    ampgradfile = [ampgrad_output_path,'/',eventid,'_ampgrad_',parameters.component,'.mat'];
    if is_eikonal_ampgrad && exist(ampgradfile,'file')
        temp = load(ampgradfile);
        ampgrad = temp.ampgrad;
    end
	if length(eventphv) ~= length(eventcs.avgphv)
		disp('Inconsist of period number for CS file and eikonal file');
		continue;
	end

	for ip = 1:length(eventphv)
		%% fit the amplitude surface
		% reset the arrays
		clear stlas stlos amps
		stlas = eventcs.stlas;
		stlos = eventcs.stlos;
		stnms = eventcs.stnms;
		if exist('badstnms','var')
			list_badstaids = find(ismember(eventcs.stnms,badstnms));
		else
			list_badstaids = [];
		end
		amps = zeros(1,length(stlas));
		for ista = 1:length(eventcs.autocor)
			if eventcs.autocor(ista).exitflag(ip)>0
				amps(ista) = eventcs.autocor(ista).amp(ip);
			else
				amps(ista) = NaN;
			end
		end
		% change from power spectrum to amplitude
		amps = amps.^.5;
        
        % Correct amplitude for local receiver effects
        if is_receiver_terms==1
            Amp_rec = receiver(ip).Amp_rec;
            for ista = 1:length(stnms)
                Istation = find(strcmp(stnms(ista),receiver(ip).stas));
                if isempty(Istation)
                    disp(['No station term for ',stnms(ista)]);
                    continue
                end
                amps(ista) = amps(ista) ./ Amp_rec(Istation);
            end
        end

		% get rid of bad stations
		badstaids = find(isnan(amps));
		stlas(badstaids) = [];
		stlos(badstaids) = [];
		amps(badstaids) = [];
        stnms(badstaids) = [];
		badstanum = 0; badstaids = [];
		for ista = 1:length(amps)
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
			meanamp = median(amps(nearstaids));
% 			if amps(ista) < meanamp./amp_var_tol | amps(ista) > meanamp.*amp_var_tol
            if amps(ista) < meanamp.*(1-amp_var_tol) | amps(ista) > meanamp.*(1+amp_var_tol)
				badstanum = badstanum+1;
				badstaids(badstanum) = ista;
			end
		end
		stlas(badstaids) = [];
		stlos(badstaids) = [];
		amps(badstaids) = [];
        stnms(badstaids) = [];
        
        if length(amps(~isnan(amps)))<min_sta_num
            ampmap = nan(size(mesh_xi));
        else
            [ampmap,mesh_xi,mesh_yi]=gridfit_jg(stlas,stlos,amps,xnode,ynode,...
                                'smooth',2,'regularizer','del4','solver','normal');
        end

		%% Calculate the correction term
        
        if is_eikonal_ampgrad == 1
            amp_grad = ampgrad(ip).dAmp';
            amp_gradlat = -ampgrad(ip).dAmpx;
            amp_gradlon = -ampgrad(ip).dAmpy;
            [~,amp_laplat,~]=delm(xi,yi,amp_gradlat);
            [~,~,amp_laplon]=delm(xi,yi,amp_gradlon);
            dAmp = amp_laplat + amp_laplon;
            dAmp = dAmp';
            amp_gradlat = amp_gradlat';
            amp_gradlon = amp_gradlon';
        else
            [amp_grad,amp_gradlat,amp_gradlon]=delm(mesh_xi,mesh_yi,ampmap);
            dAmp=del2m(mesh_xi,mesh_yi,ampmap);
        end
%         dAmp_test = del2m(mesh_xi,mesh_yi,ampmap);
%         inan = isnan(dAmp);
%         dAmp_test(inan) = nan;
%         figure(1); clf;
%         subplot(1,2,1);
%         imagesc(dAmp);
%         cb = colorbar;
%         subplot(1,2,2);
%         imagesc(dAmp_test);
%         colorbar;
%         caxis(cb.Limits);
        
        
		amp_term=-dAmp./ampmap./(2*pi/periods(ip)).^2;
		% smooth the correction term 
		smD=max([300 periods(ip).*parameters.refv]);
        if length(amps(~isnan(amps)))>min_sta_num
            amp_term = gridfit_jg(mesh_xi(:),mesh_yi(:),amp_term(:),xnode,ynode,...
                                'smooth',floor(smD./deg2km(gridsize)),'regularizer','laplacian','solver','normal');
        end
		% prepare the avg phase velocity and event phase velocity
		avgGV = avgphv(ip).GV;
		if sum(size(avgGV)==size(xi)) < 2
			avgGV = interp2(avgphv(ip).xi,avgphv(ip).yi,avgphv(ip).GV,xi,yi,'linear',NaN);
		end
		eventGV = eventphv(ip).GV;
		if sum(size(eventGV)==size(xi)) < 2
			eventxnode = eventphv(ip).lalim(1):eventphv(ip).gridsize:eventphv(ip).lalim(2);
			eventynode = eventphv(ip).lolim(1):eventphv(ip).gridsize:eventphv(ip).lolim(2);
			[eventxi eventyi] = ndgrid(eventxnode,eventynode);
			eventGV = interp2(eventxi,eventyi,eventphv(ip).GV,xi,yi,'linear',NaN);
		end
		% remove the region with small amplitude area
%		meanamp = nanmean(ampmap(:));
%		Tampmap = ampmap';
%		ind = find(Tampmap(:)<meanamp.*parameters.min_amp_tol);
%		eventGV(ind) = NaN;
		% apply correction
		[GV_cor alpha_errs alphas] = amp_correct(avgGV, eventGV, amp_term, alpha_range, alpha_search_grid);
        
        % JBR - line to replace negative phase velocities (imaginary) with event velocity
        GV_cor(imag(GV_cor)~=0) = eventGV(imag(GV_cor)~=0);
        
		[temp bestalphai] = min(alpha_errs);
		bestalpha = alphas(bestalphai);
		fprintf('%f ',bestalpha);

		% fill in informations
		helmholtz(ip).evla = eventphv(ip).evla;
		helmholtz(ip).evlo = eventphv(ip).evlo;
        helmholtz(ip).Mw = eventphv(ip).Mw;
		helmholtz(ip).raydense = eventphv(ip).raydense;
		helmholtz(ip).goodnum = eventphv(ip).goodnum;
		helmholtz(ip).badnum = eventphv(ip).badnum;
		helmholtz(ip).id = eventphv(ip).id;
		helmholtz(ip).xi = xi;
		helmholtz(ip).yi = yi;
		helmholtz(ip).GV_cor = GV_cor;
		helmholtz(ip).GV = eventGV;
		helmholtz(ip).alpha_errs = alpha_errs;
		helmholtz(ip).alphas = alphas;
		helmholtz(ip).bestalpha = bestalpha;
        helmholtz(ip).amps = amps;
		helmholtz(ip).amp_term = amp_term';
        helmholtz(ip).ampmap = ampmap';
        helmholtz(ip).amp_grad = amp_grad';
        helmholtz(ip).amp_gradlat = amp_gradlat';
        helmholtz(ip).amp_gradlon = amp_gradlon';
        helmholtz(ip).amp_lap = dAmp';
		helmholtz(ip).period = periods(ip);
        helmholtz(ip).stainfo.stlas = stlas;
        helmholtz(ip).stainfo.stlos = stlos;
		bestalphas(ip,ie) = bestalpha;

		% plot to check
		if isfigure
			figure(37)
			clf
                        set(gcf,'renderer','zbuffer');
			subplot(2,2,1)
			ax = worldmap(lalim, lolim);
			surfacem(xi,yi,eventGV);
            if ~isnan(nanmean(eventGV(:)))
                caxis([nanmean(eventGV(:))*(1-r) nanmean(eventGV(:))*(1+r)])
            end
			colorbar
			title('before cor');
			subplot(2,2,2)
			ax = worldmap(lalim, lolim);
			surfacem(xi,yi,GV_cor);
            if ~isnan(nanmean(eventGV(:)))
                caxis([nanmean(eventGV(:))*(1-r) nanmean(eventGV(:))*(1+r)])
            end
			colorbar
			title('after cor');
			nanind = find(isnan(eventGV(:)));
			ampmap = ampmap';
			ampmap(nanind) = NaN;
			amp_term = amp_term';
			amp_term(nanind) = NaN;
			subplot(2,2,3)
			ax = worldmap(lalim, lolim);
			surfacem(xi,yi,ampmap);
			title('amplitude map')
			if ~isempty(stlas)
                plotm(stlas,stlos,'v')
                la_gc = [];
                lo_gc = [];
                for ista = 1:length(stlas)
                    [la,lo]=track2('gc',eventphv(ip).evla,eventphv(ip).evlo,stlas(ista),stlos(ista));
                    la_gc = [la_gc; la; nan];
                    lo_gc = [lo_gc; lo; nan];
                end
                plotm(la_gc,lo_gc,'-k');
            end
            colormap(seiscmap)
			colorbar
			subplot(2,2,4)
			ax = worldmap(lalim, lolim);
			surfacem(xi,yi,amp_term);
			colorbar
			[temp bestalphai] = min(alpha_errs);
			title('correction term')
            drawnow;
            
			figure(38)
			clf
                        set(gcf,'renderer','zbuffer');
			plot(alphas,alpha_errs,'x');
            drawnow;
            
%             pause;
		end % end of isfigure
	end  % loop of period
    
    if is_save_amp_fig
        clim = {};
        for ip = 1:length(eventphv)
            figure(39);
            ampmap = helmholtz(ip).ampmap;
            ampmap = ampmap';
			ampmap(nanind) = NaN;
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
            la_gc = [];
            lo_gc = [];
            for ista = 1:length(helmholtz(ip).stainfo.stlas)
                [la,lo]=track2('gc',eventphv(ip).evla,eventphv(ip).evlo,helmholtz(ip).stainfo.stlas(ista),helmholtz(ip).stainfo.stlos(ista));
                la_gc = [la_gc; la; nan];
                lo_gc = [lo_gc; lo; nan];
            end
            N=3; M = floor(length(periods)/N)+1;
            subplot(M,N,ip)
            ax = worldmap(lalim, lolim);
            % Plot as percent of median
            ampmap_perc = (ampmap-nanmedian(ampmap(:)))./nanmedian(ampmap(:))*100;
%                 surfacem(xi,yi,ampmap);
            surfacem(xi,yi,ampmap_perc);
            plotm(helmholtz(ip).stainfo.stlas,helmholtz(ip).stainfo.stlos,'v');
            plotm(la_gc,lo_gc,'-k');
            title([num2str(periods(ip)),' s'],'fontsize',15)
            caxis([-40 40]);
            cb = colorbar;
            clim{ip} = cb.Limits;
            colormap(seiscmap)
        end
        for ip = 1:length(eventphv)
            figure(40);
            ampmap = helmholtz(ip).ampmap;
            ampmap = ampmap';
			ampmap(nanind) = NaN;
            amps = helmholtz(ip).amps;
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
            % Plot as percent of median
            amps_perc = (amps-nanmedian(ampmap(:)))./nanmedian(ampmap(:))*100;
%                 scatterm(stlas,stlos,100,amps,'v','filled','markeredgecolor',[0 0 0]);
            scatterm(helmholtz(ip).stainfo.stlas,helmholtz(ip).stainfo.stlos,100,amps_perc,'v','filled','markeredgecolor',[0 0 0]);
            plotm(la_gc,lo_gc,'-k');
            title([num2str(periods(ip)),' s'],'fontsize',15)
            colorbar
            caxis(clim{ip});
            colormap(seiscmap)
        end
    end
    
	matfilename = fullfile(helmholtz_path,[eventphv(1).id,'_helmholtz_',parameters.component,'.mat']);
	save(matfilename,'helmholtz');
	fprintf('\n');
	disp(['Saved to ',matfilename]);
    if is_save_amp_fig
        figdir = [workingdir,'/figs/helmholtz/'];
        if ~exist(figdir)
            mkdir(figdir);
        end
        save2pdf([figdir,eventphv(1).id,'_helmholtz_',parameters.component,'.pdf'],39,100);
        save2pdf([figdir,eventphv(1).id,'_helmholtz_',parameters.component,'_StaAmps.pdf'],40,100);
    end
end  % loop of events
