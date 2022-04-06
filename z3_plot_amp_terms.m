% Estimate travel time surface and gradient as well as amplitude gradient
% for use in eq (4) of Bao et al. (2016) GJI
%
% This version determines uncertainties via bootstrapping
%
% github.com/jbrussell
% 2021-06

clear;
% setup parameters
setup_parameters

% isoverwrite = 1;
is_figures_aux = 0; % lots of additional figures...
% is_save_amp_fig = 1;

min_Mw = parameters.min_Mw_alpha; % minimum magnitude
min_Ngrcells = parameters.min_Ngrcells; % minimum numbe of grid cells required in order to use event

is_eikonal_ampgrad_norm = parameters.is_eikonal_ampgrad_norm;

r = 0.05;

workingdir = parameters.workingdir;
eventcs_path = [workingdir,'CSmeasure/'];
phase_v_path = [workingdir,'eikonal/'];
helmholtz_path = [workingdir,'helmholtz/'];
helmholtz_stack_file = [workingdir,'helmholtz_stack_',parameters.component];
traveltime_path = [workingdir,'traveltime/'];
attenuation_path = workingdir; %[workingdir,'attenuation/'];

% if ~exist(attenuation_path,'dir')
% 	mkdir(attenuation_path);
% end

% load stacked phase velocity map
load(helmholtz_stack_file);

% set up useful variables
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
xnode = lalim(1):gridsize:lalim(2);
ynode = lolim(1):gridsize:lolim(2);
[xi yi] = ndgrid(xnode,ynode);
Nx=length(xnode); Ny=length(ynode);
alpha_range = parameters.alpha_range;
alpha_search_grid = parameters.alpha_search_grid;
periods = parameters.periods;

eventfiles = dir([traveltime_path,'/*_traveltime_',parameters.component,'.mat']);

load seiscmap

if exist('badampsta.lst','file')
	badstnms = textread('badampsta.lst','%s');
	disp('Found Bad amplitude stations:')
	for ista = 1:length(badstnms)
	disp(badstnms(ista))
	end
end

clear attenuation
for ie = 1:length(eventfiles)
    clear ampgradR_ampnorm_dot_tpgrad amp_gradR_ampnorm_map amp_gradR_map amp_gradT_map amp_gradlat_ampnorm_map amp_gradlon_ampnorm_map ampgrad_dot_tpgrad_ampnorm amp_term amp_term_err azi amp_decay_map tp_focus_map tp_grad_map amp_grad_map ampgrad_dot_tpgrad amp_grad_norm_map evids dist_map amp_map amp_gradlat_map amp_gradlon_map tp_gradlat_map tp_gradlon_map
    evcnt = 0;
    
    % read in data for this event
    clear eventphv eventcs traveltime
    load(fullfile(traveltime_path,eventfiles(ie).name));
    eventid = traveltime(1).id;
    disp(eventid);    
    eventcsfile = [eventcs_path,'/',eventid,'_cs_',parameters.component,'.mat'];
    if exist(eventcsfile,'file')
        load(eventcsfile);
    else
        disp(['Cannot find CS file for ',eventid,', Skipped']);
        continue;
    end
    eventeikonalfile = [phase_v_path,'/',eventid,'_eikonal_',parameters.component,'.mat'];
    if exist(eventeikonalfile,'file')
        load(eventeikonalfile);
    else
        disp(['Cannot find eikonal file for ',eventid,', Skipped']);
        continue;
    end
    helmholtzfile = [helmholtz_path,'/',eventid,'_helmholtz_',parameters.component,'.mat'];
    if exist(helmholtzfile,'file')
        temp = load(helmholtzfile);
        helmholtz = temp.helmholtz;
    else
        disp(['Cannot find Helmholtz file for ',eventid,', Skipped']);
        continue;
    end
    
    amp_decay_map = nan(size(xi,1),size(xi,2),length(periods));
    tp_focus_map = amp_decay_map;
    amp_term = amp_decay_map;
    for ip = 1:length(periods)
    %for ie = 59
        
        if traveltime(1).Mw < min_Mw
            continue
        end
        
        
        Inotnan = find(~isnan(traveltime(ip).GV_cor));
        Inotnan_tp = find(~isnan(traveltime(ip).tp_grad));
%         Inotnan = find(~isnan(azi_prop));
        if length(Inotnan) < min_Ngrcells || length(Inotnan_tp) < min_Ngrcells
            continue
        end
        
        evcnt = evcnt+1;

		%% fit the amplitude surface
		% reset the arrays
		clear stlas stlos tp
		stlas = eventcs.stlas;
		stlos = eventcs.stlos;
		stnms = eventcs.stnms;
        
        % Load amplitude fields
        amp = helmholtz(ip).ampmap;
        amp(isnan(traveltime(ip).GV_cor)) = nan;
        amp_grad = helmholtz(ip).amp_grad;
        amp_gradlat = helmholtz(ip).amp_gradlat;
        amp_gradlon = helmholtz(ip).amp_gradlon;
        if is_eikonal_ampgrad_norm
            amp_grad_ampnorm = helmholtz(ip).amp_grad_ampnorm;
            amp_gradlat_ampnorm = helmholtz(ip).amp_gradlat_ampnorm;
            amp_gradlon_ampnorm = helmholtz(ip).amp_gradlon_ampnorm;
            amp_grad_ampnorm_err = full(helmholtz(ip).amp_grad_ampnorm_err);
            amp_gradlat_ampnorm_err = full(helmholtz(ip).amp_gradlat_ampnorm_err);
            amp_gradlon_ampnorm_err = full(helmholtz(ip).amp_gradlon_ampnorm_err);
        end
        
        % Load travel-time fields
        tp_grad = traveltime(ip).tp_grad;
        tp_gradlat = traveltime(ip).tp_gradlat;
        tp_gradlon = traveltime(ip).tp_gradlon;
        tp_lap = traveltime(ip).tp_lap;
        tp_grad_err = full(traveltime(ip).tp_grad_err);
        tp_gradlat_err = full(traveltime(ip).tp_gradlat_err);
        tp_gradlon_err = full(traveltime(ip).tp_gradlon_err);
        tp_lap_err = full(traveltime(ip).tp_lap_err);
        
        % Get propagation azimuth at each grid cell
        azi_prop = traveltime(ip).tp_ang;
        azi_prop(azi_prop<0) = azi_prop(azi_prop<0)+360;
        
        % Get structural phase velocity
        phv = avgphv(ip).GV_cor ;
        phv(isnan(traveltime(ip).GV_cor)) = nan;
        phv_err = full(traveltime(ip).phv_err);

        % Calculate terms from Bao et al. (2016) equation 4
        % Amplitude decay term
        if is_eikonal_ampgrad_norm
            amp_decay = 2*(amp_gradlat_ampnorm.*tp_gradlat + amp_gradlon_ampnorm.*tp_gradlon);
            % Propagate errors
            amp_decay_err = 2*( (tp_gradlat.*amp_gradlat_ampnorm_err).^2 + (amp_gradlat_ampnorm.*tp_gradlat_err).^2 ...
                               +(tp_gradlon.*amp_gradlon_ampnorm_err).^2 + (amp_gradlon_ampnorm.*tp_gradlon_err).^2 ).^0.5;
        else
            amp_decay = 2*(amp_gradlat.*tp_gradlat + amp_gradlon.*tp_gradlon) ./ amp;
            amp_decay_err = nan(size(amp_decay));
        end
        % Focusing correction term
        tp_focus = tp_lap;
        tp_focus_err = tp_lap_err;
        
        % smooth the terms
        for ii=1
            if isempty(find(~isnan(amp_decay))) || isempty(find(~isnan(tp_focus))) ...
               || length(find(~isnan(amp_decay)))<5 || length(find(~isnan(tp_focus)))<5 ...
			   || length(find(~isnan(amp_decay_err)))<5 || length(find(~isnan(tp_focus_err)))<5
                continue
            end
            smD=max([300 periods(ip).*parameters.refv]);
            amp_decay = gridfit_jg(xi(:),yi(:),amp_decay(:),xnode,ynode,...
                                'smooth',floor(smD./deg2km(gridsize)),'regularizer','laplacian','solver','normal')';
            tp_focus = gridfit_jg(xi(:),yi(:),tp_focus(:),xnode,ynode,...
                                'smooth',floor(smD./deg2km(gridsize)),'regularizer','laplacian','solver','normal')';
            amp_decay_err = gridfit_jg(xi(:),yi(:),amp_decay_err(:),xnode,ynode,...
                                'smooth',floor(smD./deg2km(gridsize)),'regularizer','laplacian','solver','normal')';
            tp_focus_err = gridfit_jg(xi(:),yi(:),tp_focus_err(:),xnode,ynode,...
                                'smooth',floor(smD./deg2km(gridsize)),'regularizer','laplacian','solver','normal')';
        end
        Inan = isnan(phv) | isnan(tp_grad);
        amp_decay(Inan) = nan;
        tp_focus(Inan) = nan;
        amp_decay_err(Inan) = nan;
        tp_focus_err(Inan) = nan;

%         % Average all grid cells
%         amp_decay = nanmean(amp_decay(:)) * ones(size(amp_decay));
%         tp_focus = nanmean(tp_focus(:)) * ones(size(tp_focus));
        
        % Corrected amplitude decay
        corr_amp_decay = amp_decay + tp_focus;
        amp_decay_map(:,:,evcnt) = amp_decay;
        tp_focus_map(:,:,evcnt) = tp_focus;
        
        amp_term(:,:,evcnt) = (phv/2) .* corr_amp_decay;     
        % Propagate errors
        amp_term_err(:,:,evcnt) = ( (0.5.*(amp_decay + tp_focus).*phv_err).^2 + (0.5.*phv.*amp_decay_err).^2 + (0.5.*phv.*tp_focus_err).^2 ).^0.5;
        azi(:,:,evcnt) = azi_prop;
        
        tp_grad_map(:,:,evcnt) = tp_grad;
        amp_grad_map(:,:,evcnt) = amp_grad;
        ampgrad_dot_tpgrad(:,:,evcnt) = amp_gradlat.*tp_gradlat + amp_gradlon.*tp_gradlon;
        evids{evcnt} = eventid;
        dists = km2deg(distance(eventcs.evla,eventcs.evlo,xi,yi,referenceEllipsoid('GRS80'))/1000);
        dist_map(:,:,evcnt) = dists;
        amp_map(:,:,evcnt) = amp;
        amp_gradlat_map(:,:,evcnt) = amp_gradlat;
        amp_gradlon_map(:,:,evcnt) = amp_gradlon;
        amp_gradR_map(:,:,evcnt)  = amp_gradlat.*cosd(azi_prop) + amp_gradlon.*sind(azi_prop);
        amp_gradT_map(:,:,evcnt)  = -amp_gradlat.*sind(azi_prop) + amp_gradlon.*cosd(azi_prop);
        tp_gradlat_map(:,:,evcnt) = tp_gradlat;
        tp_gradlon_map(:,:,evcnt) = tp_gradlon;
        if is_eikonal_ampgrad_norm
            amp_gradlat_ampnorm_map(:,:,evcnt) = amp_gradlat_ampnorm;
            amp_gradlon_ampnorm_map(:,:,evcnt) = amp_gradlon_ampnorm;
            amp_grad_norm_map(:,:,evcnt) = amp_grad_ampnorm;
            ampgrad_dot_tpgrad_ampnorm(:,:,evcnt) = amp_gradlat_ampnorm.*tp_gradlat + amp_gradlon_ampnorm.*tp_gradlon;
            amp_gradR_ampnorm_map(:,:,evcnt) = amp_gradlat_ampnorm.*cosd(azi_prop) + amp_gradlon_ampnorm.*sind(azi_prop);
            ampgradR_ampnorm_dot_tpgrad(:,:,evcnt) = amp_gradR_ampnorm_map(:,:,evcnt) .* tp_grad;
        else
            amp_gradlat_ampnorm_map(:,:,evcnt) = amp_gradlat ./ amp;
            amp_gradlon_ampnorm_map(:,:,evcnt) = amp_gradlon ./ amp;
            amp_grad_norm_map(:,:,evcnt) = amp_grad ./ amp;
            ampgrad_dot_tpgrad_ampnorm(:,:,evcnt) = ampgrad_dot_tpgrad(:,:,evcnt) ./ amp;
            amp_gradR_ampnorm_map(:,:,evcnt) = amp_gradR_map(:,:,evcnt) ./ amp;
            ampgradR_ampnorm_dot_tpgrad(:,:,evcnt) = (amp_gradR_map(:,:,evcnt) ./ amp) .* tp_grad;
        end
    end
    
    figure(40); clf;
    set(gcf,'Position',[84           3         744        1022]);
    sgtitle(['$\frac{2\nabla A \cdot \nabla \tau}{A}$  ',eventphv(ip).id,' M',num2str(eventphv(ip).Mw)],'fontweight','bold','fontsize',18,'interpreter','latex')
    N=3; M = floor(length(periods)/N)+1;
    clim = {};
    for ip = 1:length(periods)    
        nanind = find(isnan(helmholtz(ip).GV(:)));
        amp_decay_map_pl = amp_decay_map(:,:,ip);
        amp_decay_map_pl(nanind) = NaN;
        [la_gc,lo_gc]=track2('gc',helmholtz(ip).evla,helmholtz(ip).evlo,mean(lalim),mean(lolim));
        subplot(M,N,ip); hold on;
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,amp_decay_map_pl);
        plotm(la_gc,lo_gc,'-k');
        title([num2str(periods(ip)),' s'],'fontsize',15)
        if ~isempty(amp_decay_map_pl(~isnan(amp_decay_map_pl(:))))
            caxis([-1 1]*max(abs(amp_decay_map_pl(:))));
        end
        cb = colorbar;
        colormap(seiscmap)
        clim{ip} = cb.Limits;
    end
    
    figure(41); clf;
    set(gcf,'Position',[84           3         744        1022]);
    sgtitle(['$\nabla^2 \tau$  ',eventphv(ip).id,' M',num2str(eventphv(ip).Mw)],'fontweight','bold','fontsize',18,'interpreter','latex')
    N=3; M = floor(length(periods)/N)+1;
    for ip = 1:length(periods)    
        nanind = find(isnan(helmholtz(ip).GV(:)));
        tp_focus_map_pl = tp_focus_map(:,:,ip);
        tp_focus_map_pl(nanind) = NaN;
        [la_gc,lo_gc]=track2('gc',helmholtz(ip).evla,helmholtz(ip).evlo,mean(lalim),mean(lolim));
        subplot(M,N,ip); hold on;
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,tp_focus_map_pl);
        plotm(la_gc,lo_gc,'-k');
        title([num2str(periods(ip)),' s'],'fontsize',15)
        colorbar
        colormap(seiscmap)
        caxis(clim{ip});
    end
    
    figure(42); clf;
    set(gcf,'Position',[84           3         744        1022]);
    sgtitle(['Corrected Amplitude Decay  ',eventphv(ip).id,' M',num2str(eventphv(ip).Mw)],'fontweight','bold','fontsize',18,'interpreter','latex')
    N=3; M = floor(length(periods)/N)+1;
    for ip = 1:length(periods)    
        nanind = find(isnan(helmholtz(ip).GV(:)));
%         amp_term_pl = amp_term(:,:,ip);
        amp_term_pl = amp_decay_map(:,:,ip) + tp_focus_map(:,:,ip);
        amp_term_pl(nanind) = NaN;
        [la_gc,lo_gc]=track2('gc',helmholtz(ip).evla,helmholtz(ip).evlo,mean(lalim),mean(lolim));
        subplot(M,N,ip); hold on;
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,amp_term_pl);
        plotm(la_gc,lo_gc,'-k');
        title([num2str(periods(ip)),' s'],'fontsize',15)
        colorbar
        colormap(seiscmap)
        caxis(clim{ip});
    end
    
    figdir = [workingdir,'/figs/ampterms/'];
    if ~exist(figdir)
        mkdir(figdir);
    end
    save2pdf([figdir,helmholtz(1).id,'_ampterms_',parameters.component,'_ApparAmpDecay.pdf'],40,100);
    save2pdf([figdir,helmholtz(1).id,'_ampterms_',parameters.component,'_FocusingCorr.pdf'],41,100);
%     save2pdf([figdir,helmholtz(1).id,'_ampterms_',parameters.component,'_ApparAlpha.pdf'],42,100);
    save2pdf([figdir,helmholtz(1).id,'_ampterms_',parameters.component,'_CorrAmpDecay.pdf'],42,100);
end

