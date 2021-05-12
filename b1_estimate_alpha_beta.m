% Estimate travel time surface and gradient as well as amplitude gradient
% for use in eq (4) of Bao et al. (2016) GJI
% github.com/jbrussell
% 2021-05

clear;

isoverwrite = 0;
isfigure = 1;
is_save_amp_fig = 1;

min_Mw = 5.5; % minimum magnitude
min_Ngrcells = 20; % minimum numbe of grid cells required in order to use event
azi_bin_deg = 30; % [deg] size of azimuthal bins
min_nbin = 10; % minimum number of measurements in order to include bin
N_min_evts = 10; % minimum number of events contributing to grid cell required in order to be considered

% setup parameters
setup_parameters

r = 0.05;

workingdir = parameters.workingdir;
eventcs_path = [workingdir,'CSmeasure/'];
phase_v_path = [workingdir,'eikonal/'];
helmholtz_path = [workingdir,'helmholtz/'];
helmholtz_stack_file = [workingdir,'helmholtz_stack_',parameters.component];
traveltime_path = [workingdir,'traveltime/'];

if ~exist(traveltime_path,'dir')
	mkdir(traveltime_path);
end

% load stacked phase velocity map
load(helmholtz_stack_file);

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
for ip = 1:length(avgphv)
    clear amp_term azi
    evcnt = 0;
    for ie = 1:length(eventfiles)
    %for ie = 59
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
        
        if traveltime(1).Mw < min_Mw
            continue
        end
        
%         % Remove grid cells with outlier propagation azimuth
%         max_degrelmean = 20;
%         Inan = isnan(traveltime(ip).GV_cor);
%         azi_prop = 90 - atan2d(traveltime(ip).tp_gradlat',traveltime(ip).tp_gradlon');
%         azi_prop(azi_prop<0) = azi_prop(azi_prop<0)+360;
%         azi_prop(Inan) = nan;
% %         Ibad_prop = abs(azi_prop-nanmean(azi_prop(:))) > max_degrelmean;
% %         azi_prop(Ibad_prop) = nan;      
%         azi_GV = angle(eventphv(ip).GVx + eventphv(ip).GVy.*sqrt(-1));
% 		azi_GV = rad2deg(azi_GV) + 180;
%         absDiffDeg = @(a,b) abs(diff(unwrap([a,b]/180*pi)*180/pi));
%         for ii=1:length(xnode)
%             for jj=1:length(ynode)
%                 diffdeg(ii,jj) = absDiffDeg(azi_GV(ii,jj),azi_prop(ii,jj));
%             end
%         end
%         Ibad_prop = diffdeg > max_degrelmean;
%         avgphv(ip).GV_cor(Ibad_prop) = nan;
%         traveltime(ip).GV_cor(Ibad_prop) = nan;
        
        Inotnan = find(~isnan(traveltime(ip).GV_cor));
%         Inotnan = find(~isnan(azi_prop));
        if length(Inotnan) < min_Ngrcells
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
        amp_grad = helmholtz(ip).amp_grad;
        amp_gradlat = helmholtz(ip).amp_gradlat;
        amp_gradlon = helmholtz(ip).amp_gradlon;
        
        % Load travel-time fields
        tp_grad = traveltime(ip).tp_grad;
        tp_gradlat = traveltime(ip).tp_gradlat;
        tp_gradlon = traveltime(ip).tp_gradlon;
        tp_lap = traveltime(ip).tp_lap;
        
        % Get propagation azimuth at each grid cell
        azi_prop = traveltime(ip).tp_ang;
        azi_prop(azi_prop<0) = azi_prop(azi_prop<0)+360;
        
        % Get structural phase velocity
        phv = avgphv(ip).GV_cor ;
        phv(isnan(traveltime(ip).GV_cor)) = nan;

%         amp(Inan)=nan; amp_grad(Inan)=nan; amp_gradlat(Inan)=nan; amp_gradlon(Inan)=nan;
%         tp_lap(Inan)=nan; tp_grad(Inan)=nan; tp_gradlat(Inan)=nan; tp_gradlon(Inan)=nan;
        
%         figure(1);
%         subplot(1,2,1);
%         tp_grad(isnan(traveltime(ip).GV)) = nan;
%         imagesc(tp_grad);
%         cb = colorbar;
%         subplot(1,2,2);
%         imagesc(1./traveltime(ip).GV);
%         colorbar;
%         caxis(cb.Limits);

        % Calculate terms from Bao et al. (2016) equation 4
        % Amplitude decay term
        amp_decay = 2*(amp_gradlat.*tp_gradlat + amp_gradlon.*tp_gradlon) ./ amp;
        % Focusing correction term
        tp_focus = tp_lap;
        
        % smooth the terms
        for ii=1
			if isempty(find(~isnan(amp_decay))) || isempty(find(~isnan(tp_focus)))
				continue
			end
            smD=max([300 periods(ip).*parameters.refv]);
            amp_decay = gridfit_jg(xi(:),yi(:),amp_decay(:),xnode,ynode,...
                                'smooth',floor(smD./deg2km(gridsize)),'regularizer','laplacian','solver','normal')';
            tp_focus = gridfit_jg(xi(:),yi(:),tp_focus(:),xnode,ynode,...
                                'smooth',floor(smD./deg2km(gridsize)),'regularizer','laplacian','solver','normal')';
        end
        Inan = isnan(phv);
        amp_decay(Inan) = nan;
        tp_focus(Inan) = nan;
        
        % Corrected amplitude decay
        corr_amp_decay = amp_decay + tp_focus;
        
        amp_term(:,:,evcnt) = (phv/2) .* corr_amp_decay;        
        azi(:,:,evcnt) = azi_prop;
    end

%     % Select central grid points
%     nanmat = nan(length(xnode),length(ynode),evcnt);
% %     nanmat(9,9,:) = 1;
%     nanmat([9-1:9+1],[9-1:9+1],:) = 1;
%     amp_term = amp_term .* nanmat;
%     azi = azi .* nanmat;
    
    
    % Remove large outliers
    Ibad = amp_term>abs(nanmedian(amp_term(:)))*100 | amp_term<abs(nanmedian(amp_term(:)))*-100;
    amp_term(Ibad) = nan;
    
    % Bin measurements by azimuth
    [~,Isort] = sort(azi(:));
    azi_srt = azi(Isort);
    amp_term_srt = amp_term(Isort);
    bins = [0:azi_bin_deg:360];
    amp_bin = 0; azi_bin = 0; amp_bin_err = 0; azi_bin_err = 0;
    for ibin = 1:length(bins)-1
        I_bin = azi_srt>=bins(ibin) & azi_srt<bins(ibin+1);
        if sum(I_bin)<min_nbin
            I_bin = false(size(I_bin));
        end
        % Weighted means and standard deviations
        amp_bin(ibin) = nanmedian(amp_term_srt(I_bin));
        amp_bin_std(ibin) = nanstd(amp_term_srt(I_bin));
        azi_bin(ibin) = (bins(ibin)+bins(ibin+1))/2;
    end
    
    %% Do curve fitting (eq 9 in Bao et al. 2016)
    
    % Unbinned 1-D sinusoidal fit
    [para, alpha, dlnbeta_dx, dlnbeta_dy]=fit_alpha_beta(azi(:),amp_term(:));
    parastd=confint(para,.95);
    dlnbeta_dx_err = (parastd(2,1)-parastd(1,1))/2;
    dlnbeta_dy_err = (parastd(2,2)-parastd(1,2))/2;
    alpha_err = (parastd(2,3)-parastd(1,3))/2;
    attenuation(ip).alpha_1d = alpha;
    attenuation(ip).dlnbeta_dx_1d = dlnbeta_dx;
    attenuation(ip).dlnbeta_dy_1d = dlnbeta_dy;
    attenuation(ip).alpha_1d_err = alpha_err;
    attenuation(ip).dlnbeta_dx_1d_err = dlnbeta_dx_err;
    attenuation(ip).dlnbeta_dy_1d_err = dlnbeta_dy_err;
    attenuation(ip).para_1d = para;
    
    % Binned 1-D sinusoidal fit
    [para, alpha, dlnbeta_dx, dlnbeta_dy]=fit_alpha_beta(azi_bin(:),amp_bin(:));
    parastd=confint(para,.95);
    dlnbeta_dx_err = (parastd(2,1)-parastd(1,1))/2;
    dlnbeta_dy_err = (parastd(2,2)-parastd(1,2))/2;
    alpha_err = (parastd(2,3)-parastd(1,3))/2;
    attenuation(ip).alpha_1d_bin = alpha;
    attenuation(ip).dlnbeta_dx_1d_bin = dlnbeta_dx;
    attenuation(ip).dlnbeta_dy_1d_bin = dlnbeta_dy;
    attenuation(ip).alpha_1d_bin_err = alpha_err;
    attenuation(ip).dlnbeta_dx_1d_bin_err = dlnbeta_dx_err;
    attenuation(ip).dlnbeta_dy_1d_bin_err = dlnbeta_dy_err;
    attenuation(ip).para_1d_bin = para;
    
    % Unbinned 1-D azimuthal average
    alpha_1d_avg = nanmean(-amp_term(:));
    attenuation(ip).alpha_1d_avg = alpha_1d_avg;
    attenuation(ip).alpha_1d_avg_err = nanstd(-amp_term(:)-alpha_1d_avg);
    
    % Unbinned 2-D sinusoidal fit
    for ix = 1:length(xnode)
        for iy = 1:length(ynode)
            amps = squeeze(amp_term(ix,iy,:));
            azis = squeeze(azi(ix,iy,:));
            if length(find(~isnan(amps)))>N_min_evts
                [para, alpha, dlnbeta_dx, dlnbeta_dy]=fit_alpha_beta(azis,amps);
                parastd=confint(para,.95);
                dlnbeta_dx_err = (parastd(2,1)-parastd(1,1))/2;
                dlnbeta_dy_err = (parastd(2,2)-parastd(1,2))/2;
                alpha_err = (parastd(2,3)-parastd(1,3))/2;
            else
                alpha = nan;
                dlnbeta_dx = nan;
                dlnbeta_dy = nan;
                alpha_err = nan;
                dlnbeta_dx_err = nan;
                dlnbeta_dy_err = nan;
                para = [];
            end            
            attenuation(ip).alpha_2d(ix,iy) = alpha;
            attenuation(ip).dlnbeta_dx_2d(ix,iy) = dlnbeta_dx;
            attenuation(ip).dlnbeta_dy_2d(ix,iy) = dlnbeta_dy;
            attenuation(ip).alpha_2d_err(ix,iy) = alpha_err;
            attenuation(ip).dlnbeta_dx_2d_err(ix,iy) = dlnbeta_dx_err;
            attenuation(ip).dlnbeta_dy_2d_err(ix,iy) = dlnbeta_dy_err;
            attenuation(ip).para_2d{ix,iy} = para;
            attenuation(ip).amp_term_2d = amp_term;
            attenuation(ip).azi = azi;
            attenuation(ip).period = periods(ip);
            attenuation(ip).evcnt = evcnt;
        end
    end
    
    figure(41);
    if ip==1
        clf;
    end
    N=3; M = floor(length(periods)/N)+1;
    subplot(M,N,ip)
    plot(azi(:),amp_term(:),'o'); hold on;
    plot(azi(Ibad),amp_term(Ibad),'o');
    plot(azi_bin(:),amp_bin(:),'o');
    x = [0:360];
    pred = attenuation(ip).dlnbeta_dx_1d*sind(x)+attenuation(ip).dlnbeta_dy_1d*cosd(x)-attenuation(ip).alpha_1d;
    plot(x,pred,'-r');
    title([num2str(attenuation(ip).period),' s'])
    xlabel('azimuth');
    ylabel('amplitude term');
end

%% Plot 1D average alpha

figure(42); clf; set(gcf,'color','w');
alpha_zhitu = [4.1 7.3 8.2 8.9 6.9]*1e-5;
f_mhz_zhitu = [10 15 20 25 30];
alphas = [attenuation(:).alpha_1d];
alphas_err = [attenuation(:).alpha_1d_err];
alphas_bin = [attenuation(:).alpha_1d_bin];
alphas_bin_err = [attenuation(:).alpha_1d_bin_err];
alphas_avg = [attenuation(:).alpha_1d_avg];
alphas_avg_err = [attenuation(:).alpha_1d_avg_err];
for ip = 1:length(attenuation)
    alphas_2d(ip) = nanmean(attenuation(ip).alpha_2d(:));
%     alphas_2d_err(ip) = nanmean(attenuation(ip).alpha_2d_err(:));
    alphas_2d_err(ip) = nanstd(attenuation(ip).alpha_2d(:));
end
errorbar(periods,alphas_bin,alphas_bin_err,'-om'); hold on;
errorbar(periods,alphas,alphas_err,'-ok');
plot(periods,alphas_avg,'-oc');
% errorbar(periods,alphas_2d,alphas_2d_err,'-ob');
plot(periods,alphas_2d,'-ob');
plot(1./f_mhz_zhitu*1000,alpha_zhitu,'xr','linewidth',3,'MarkerSize',8); hold on;
plot(1./f_mhz_zhitu*1000,alpha_zhitu*2,'x','color',[0 0.85 0],'linewidth',3,'MarkerSize',8); hold on;
legend({'1D fit (bin)','1D fit','1D mean','2D fit avg','Zhitu','Zhitu x 2'},'location','northeastoutside','fontsize',15)
set(gca,'fontsize',15,'linewidth',1.5);
xlabel('Period (s)');
ylabel('\alpha (km^{-1})');
xlim([min(periods)-10 max(periods)+10]);

%% Plot 2D maps of alpha
figure(43); clf; set(gcf,'position',[146           1         726        1024],'color','w');
N=3; M = floor(length(periods)/N)+1;
for ip = 1:length(attenuation)    
    alpha_2d = attenuation(ip).alpha_2d;    
    subplot(M,N,ip)
    ax = worldmap(lalim, lolim);
    surfacem(xi,yi,alpha_2d);
%     if ~isempty(stlas) 
%         plotm(stlas,stlos,'v');
%     end
    title([num2str(periods(ip)),' s'],'fontsize',15)
    cb = colorbar;
    caxis([-1e-4 4e-4]);
    colormap(flip(seiscmap))
end


%% Save
% matfilename = fullfile(traveltime_path,[eventphv(1).id,'_attenuation_',parameters.component,'.mat']);
% save(matfilename,'attenuation');
% fprintf('\n');
% disp(['Saved to ',matfilename]);
% if is_save_amp_fig
%     figdir = [workingdir,'/figs/attenuation/'];
%     if ~exist(figdir)
%         mkdir(figdir);
%     end
%     save2pdf([figdir,eventphv(1).id,'_attenuation_',parameters.component,'.pdf'],39,100);
%     save2pdf([figdir,eventphv(1).id,'_attenuation_',parameters.component,'_StaAmps.pdf'],40,100);
% end