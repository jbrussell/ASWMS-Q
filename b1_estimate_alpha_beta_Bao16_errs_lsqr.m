% Estimate travel time surface and gradient as well as amplitude gradient
% for use in eq (4) of Bao et al. (2016) GJI
% github.com/jbrussell
% 2021-06

clear;
% setup parameters
setup_parameters

% isoverwrite = 1;
is_figures_aux = 0; % lots of additional figures...
% is_save_amp_fig = 1;
is_save_mat = 1;

% min_Mw = 5.5; % minimum magnitude
% min_Ngrcells = 20; % minimum numbe of grid cells required in order to use event
% azi_bin_deg = 30; % [deg] size of azimuthal bins
% min_nbin = 10; % minimum number of measurements in order to include bin
% N_min_evts = 10; % minimum number of events contributing to grid cell required in order to be considered

min_Mw = parameters.min_Mw_alpha; % minimum magnitude
min_Ngrcells = parameters.min_Ngrcells; % minimum numbe of grid cells required in order to use event
azi_bin_deg = parameters.azi_bin_deg; % [deg] size of azimuthal bins
min_nbin = parameters.min_nbin; % minimum number of measurements in order to include bin
N_min_evts = parameters.N_min_evts; % minimum number of events contributing to grid cell required in order to be considered
smsize_alpha = parameters.smsize_alpha; % number of nearby gridcells to gather data from
smweight_beta = parameters.smweight_beta; % Second derivative smoothing weight for beta map
smooth_alpha_Nwl = parameters.smooth_alpha_Nwl; % [wavelengths] smoothing radius of 2d alpha map

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
for ip = 1:length(avgphv)
    clear ampgradR_ampnorm_dot_tpgrad amp_gradR_ampnorm_map amp_gradR_map amp_gradT_map amp_gradlat_ampnorm_map amp_gradlon_ampnorm_map ampgrad_dot_tpgrad_ampnorm amp_term amp_term_err azi amp_decay_map tp_focus_map tp_grad_map amp_grad_map ampgrad_dot_tpgrad amp_grad_norm_map evids dist_map amp_map amp_gradlat_map amp_gradlon_map tp_gradlat_map tp_gradlon_map
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
        
%         % Test setting laplacian(tp)=0 and grad(tp)=const;
%         tp_lap = zeros(size(tp_lap));
%         tp_gradlat = nanmean(tp_gradlat(:)) * ones(size(tp_gradlat));
%         tp_gradlon = nanmean(tp_gradlon(:)) * ones(size(tp_gradlon));
%         tp_grad = sqrt(tp_gradlat.^2 + tp_gradlon.^2);
%         azi_prop = 90 - atan2d(tp_gradlat,tp_gradlon);
%         azi_prop(azi_prop<0) = azi_prop(azi_prop<0)+360;

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
	amp_term(isnan(azi(:))) = nan;
    
	% Ensure that no individual events or pixels dominate the weighting scheme
	cutoff = prctile(amp_term_err(:),5);
	amp_term_err(amp_term_err<cutoff) = cutoff;
	amp_term_err(isnan(amp_term_err)) = inf;

	% Do initial (unweighted) fitting and remove outliers
	[~, alpha, dlnbeta_dx, dlnbeta_dy]=fit_alpha_beta(azi(:),amp_term(:));
	pre = dlnbeta_dx*sind(azi) + dlnbeta_dy*cosd(azi) - alpha;
	res = amp_term(:) - pre(:);
	stdres = nanstd(res);
	ibad = find(abs(res) > 2*stdres);
	amp_term_err(ibad) = inf;
	amp_term(ibad) = nan;

%     % Select central grid points
%     nanmat = nan(length(xnode),length(ynode),evcnt);
% %     nanmat(9,9,:) = 1;
%     nanmat([9-1:9+1],[9-1:9+1],:) = 1;
%     amp_term = amp_term .* nanmat;
%     azi = azi .* nanmat;
    
    
%     % Remove large outliers
%     Ibad = amp_term>abs(nanmedian(amp_term(:)))*100 | amp_term<abs(nanmedian(amp_term(:)))*-100;
%     amp_term(Ibad) = nan;
    
    % Bin measurements by azimuth
    [~,Isort] = sort(azi(:));
    azi_srt = azi(Isort);
    amp_term_srt = amp_term(Isort);
    bins = [0:azi_bin_deg:360];
    amp_bin = 0; azi_bin = 0; amp_bin_std = 0; azi_bin_err = 0;
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
	isbad = isnan(amp_bin_std) | amp_bin_std==0;
    amp_bin(isbad) = [];
    amp_bin_std(isbad) = [];
    azi_bin(isbad) = [];
    
    %% Do curve fitting (eq 9 in Bao et al. 2016)
    attenuation(ip).evids = evids;
    attenuation(ip).amp_term = amp_term;
    attenuation(ip).amp_term_err = amp_term_err;
    attenuation(ip).azi = azi;
    attenuation(ip).dist_map = dist_map;
    attenuation(ip).amp_decay_map = amp_decay_map;
    attenuation(ip).tp_focus_map = tp_focus_map;
    attenuation(ip).tp_grad_map = tp_grad_map;
    attenuation(ip).amp_grad_map = amp_grad_map;
    attenuation(ip).amp_grad_norm_map = amp_grad_norm_map;
    attenuation(ip).ampgrad_dot_tpgrad = ampgrad_dot_tpgrad;
    attenuation(ip).amp_map = amp_map;
    attenuation(ip).amp_gradlat_map = amp_gradlat_map;
    attenuation(ip).amp_gradlon_map = amp_gradlon_map;
    attenuation(ip).tp_gradlat_map = tp_gradlat_map;
    attenuation(ip).tp_gradlon_map = tp_gradlon_map;
    attenuation(ip).amp_gradR_map = amp_gradR_map;
    attenuation(ip).amp_gradT_map = amp_gradT_map;
    attenuation(ip).amp_gradlat_ampnorm_map = amp_gradlat_ampnorm_map;
    attenuation(ip).amp_gradlon_ampnorm_map = amp_gradlon_ampnorm_map;
    attenuation(ip).ampgrad_dot_tpgrad_ampnorm = ampgrad_dot_tpgrad_ampnorm;
    attenuation(ip).amp_gradR_ampnorm_map = amp_gradR_ampnorm_map;
    attenuation(ip).ampgradR_ampnorm_dot_tpgrad = ampgradR_ampnorm_dot_tpgrad;
    
    % Unbinned 1-D sinusoidal fit
%     [para, alpha, dlnbeta_dx, dlnbeta_dy]=fit_alpha_beta(azi(:),amp_term(:),amp_term_err(:));
    [alpha, dlnbeta_dx, dlnbeta_dy, beta_tau, azi_maxamp, alpha_err, dlnbeta_dx_err, dlnbeta_dy_err, beta_tau_err, azi_maxamp_err]=fit_alpha_beta_lsqr(azi(:),amp_term(:),amp_term_err(:));
    attenuation(ip).alpha_1d = alpha;
    attenuation(ip).beta_tau_1d = beta_tau;
    attenuation(ip).azi_maxamp_1d = azi_maxamp;
    attenuation(ip).alpha_1d_err = alpha_err;
    attenuation(ip).beta_tau_1d_err = beta_tau_err;
    attenuation(ip).azi_maxamp_1d_err = azi_maxamp_err;
    
    % Binned 1-D sinusoidal fit
    [alpha, dlnbeta_dx, dlnbeta_dy, beta_tau, azi_maxamp, alpha_err, dlnbeta_dx_err, dlnbeta_dy_err, beta_tau_err, azi_maxamp_err]=fit_alpha_beta_lsqr(azi_bin(:),amp_bin(:),amp_bin_std(:));
    attenuation(ip).alpha_1d_bin = alpha;
    attenuation(ip).beta_tau_1d_bin = beta_tau;
    attenuation(ip).azi_maxamp_1d_bin = azi_maxamp;
    attenuation(ip).alpha_1d_bin_err = alpha_err;
    attenuation(ip).beta_tau_1d_bin_err = beta_tau_err;
    attenuation(ip).azi_maxamp_1d_bin_err = azi_maxamp_err;
	attenuation(ip).amp_bin = amp_bin;
	attenuation(ip).amp_bin_std = amp_bin_std;
	attenuation(ip).azi_bin = azi_bin;
    
    % Unbinned 1-D azimuthal average
    alpha_1d_avg = nanmean(-amp_term(:));
    attenuation(ip).alpha_1d_avg = alpha_1d_avg;
    attenuation(ip).alpha_1d_avg_err = nanstd(-amp_term(:)-alpha_1d_avg);
    
    % Unbinned 2-D sinusoidal fit
    for ix = 1:length(xnode)
        for iy = 1:length(ynode)
            lowx=max(1,ix-smsize_alpha);
            upx=min(Nx,ix+smsize_alpha);
            lowy=max(1,iy-smsize_alpha);
            upy=min(Ny,iy+smsize_alpha);
            amps=[]; azis=[]; amps_err=[];
            for ii=lowx:upx
                for jj=lowy:upy
                    amps = [amps; squeeze(amp_term(ii,jj,:))];
                    amps_err = [amps_err; squeeze(amp_term_err(ii,jj,:))];
                    azis = [azis; squeeze(azi(ii,jj,:))];
                end
            end
			isbin = 0;
            if isbin
				% Bin measurements by azimuth
	            [~,Isort] = sort(azis(:));
	            azi_srt = azis(Isort);
	            amp_term_srt = amps(Isort);
	            bins = [0:azi_bin_deg:360];
	            amp_bin = 0; azi_bin = 0; amp_bin_std = 0; azi_bin_err = 0;
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
                amps = amp_bin;
                azis = azi_bin;
            end
%             amps = squeeze(amp_term(ix,iy,:));
%             azis = squeeze(azi(ix,iy,:));
            if length(find(~isnan(amps)))>N_min_evts
                [alpha, dlnbeta_dx, dlnbeta_dy, beta_tau, azi_maxamp, alpha_err, dlnbeta_dx_err, dlnbeta_dy_err, beta_tau_err, azi_maxamp_err]=fit_alpha_beta_lsqr(azis(:),amps(:),amps_err(:));
            else
                alpha = nan;
                beta_tau = nan;
                dlnbeta_dx = nan;
                dlnbeta_dy = nan;
                azi_maxamp = nan;
                alpha_err = nan;
                beta_tau_err = nan;
                azi_maxamp_err = nan;
            end            
            attenuation(ip).alpha_2d(ix,iy) = alpha;
            attenuation(ip).dlnbeta_dx_2d(ix,iy) = dlnbeta_dx;
            attenuation(ip).dlnbeta_dy_2d(ix,iy) = dlnbeta_dy;
            attenuation(ip).beta_tau_2d(ix,iy) = beta_tau;
            attenuation(ip).azi_maxamp_2d(ix,iy) = azi_maxamp;
            attenuation(ip).alpha_2d_err(ix,iy) = alpha_err;
			attenuation(ip).dlnbeta_dx_2d_err(ix,iy) = dlnbeta_dx_err;
            attenuation(ip).dlnbeta_dy_2d_err(ix,iy) = dlnbeta_dy_err;
            attenuation(ip).beta_tau_2d_err(ix,iy) = beta_tau_err;
            attenuation(ip).azi_maxamp_2d_err(ix,iy) = azi_maxamp_err;
            attenuation(ip).amp_term_2d = amp_term;
            attenuation(ip).azi = azi;
            attenuation(ip).period = periods(ip);
            attenuation(ip).evcnt = evcnt;
        end
    end
	% Smooth alpha map based on wavelength
    D = smooth_alpha_Nwl*nanmean(avgphv(ip).GV_cor(:))*periods(ip);
    alpha_2d_sm = smoothmap(xi,yi,attenuation(ip).alpha_2d,D);
    alpha_2d_sm(find(isnan(attenuation(ip).alpha_2d))) = NaN;
    attenuation(ip).alpha_2d = alpha_2d_sm;
	
    % Invert for beta map
	[beta_2d,chi2] = inv_beta(xi,yi,attenuation(ip).dlnbeta_dy_2d,attenuation(ip).dlnbeta_dx_2d,attenuation(ip).dlnbeta_dy_2d_err,attenuation(ip).dlnbeta_dx_2d_err,smweight_beta);
    beta_2d(isnan(attenuation(ip).dlnbeta_dy_2d)) = nan;
    attenuation(ip).beta_2d = beta_2d;
	attenuation(ip).beta_2d_chi2 = chi2;
	
	% figure(1000); clf;
    % subplot(1,2,1); hold on;
    % ax10001 = worldmap(lalim, lolim);
    % surfacem(xi,yi,attenuation(ip).beta_tau_2d);
    % quiverm(xi,yi,attenuation(ip).dlnbeta_dy_2d,attenuation(ip).dlnbeta_dx_2d,'-k')
    % colorbar;
    % subplot(1,2,2); hold on;
    % ax10002 = worldmap(lalim, lolim);
    % surfacem(xi,yi,beta_2d);
    % quiverm(xi,yi,attenuation(ip).dlnbeta_dy_2d,attenuation(ip).dlnbeta_dx_2d,'-k')
    % colorbar;
    
    figure(41);
    if ip==1
        clf;
        set(gcf,'Position',[616    71   850   947]);
    end
    N=3; M = floor(length(periods)/N)+1;
    subplot(M,N,ip)
%     plot(azi(:),amp_term(:),'.b'); hold on;
    errorbar(azi(:),amp_term(:),amp_term_err(:),'.b'); hold on;
%     plot(squeeze(azi(9,9,:)),squeeze(amp_term(9,9,:)),'.r'); hold on;
%     plot(azi(Ibad),amp_term(Ibad),'o');
	errorbar(attenuation(ip).azi_bin(:),attenuation(ip).amp_bin(:),attenuation(ip).amp_bin_std(:),'og');
    x = [0:360];
    pred = attenuation(ip).beta_tau_1d.*cosd(x-attenuation(ip).azi_maxamp_1d)-attenuation(ip).alpha_1d;
    plot(x,pred,'-r');
    title([num2str(attenuation(ip).period),' s'])
    xlabel('azimuth');
    ylabel('amplitude term');
    ylim([-5e-4 2e-4]);
end

%% Plot 1D average alpha

mode = readMINEOS_qfile('./qfiles/pa5_5km.s0to66.q',0);
% mode = readMINEOS_qfile('./qfiles/S362ANI_NoMelt.s0to100.q',0);
alpha_MINEOS = mode.wrad ./ (2*mode.grv) ./ mode.q;

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
% Report which ever is largest, std of map or mean of stds
alphas_2d_err(ip) = max([nanstd(attenuation(ip).alpha_2d(:)), nanmean(attenuation(ip).alpha_2d_err(:))]);
end
plot(mode.T,alpha_MINEOS,'-','color',[0.7 0.7 0.7],'linewidth',5); hold on;
errorbar(periods,alphas_bin,alphas_bin_err,'-om'); hold on;
errorbar(periods,alphas,alphas_err,'-ok');   
plot(periods,alphas_avg,'-oc');
% errorbar(periods,alphas_2d,alphas_2d_err,'-ob');
plot(periods,alphas_2d,'-ob');
plot(1./f_mhz_zhitu*1000,alpha_zhitu,'xr','linewidth',3,'MarkerSize',8); hold on;
plot(1./f_mhz_zhitu*1000,alpha_zhitu*2,'x','color',[0 0.85 0],'linewidth',3,'MarkerSize',8); hold on;
legend({'True','1D fit (bin)','1D fit','1D mean','2D fit avg','Zhitu','Zhitu x 2'},'location','northeastoutside','fontsize',15)
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
    ylabel(cb,'\alpha');
    caxis([-1e-4 4e-4]);
    colormap(flip(seiscmap))
end

%% Plot 2D maps of Beta Gradient
figure(63); clf; set(gcf,'position',[146           1         726        1024],'color','w');
N=3; M = floor(length(periods)/N)+1;
for ip = 1:length(attenuation)    
    alpha_2d = attenuation(ip).alpha_2d;   
    beta_tau_2d = attenuation(ip).beta_tau_2d;
    dlnbeta_dx_2d = attenuation(ip).dlnbeta_dx_2d;
    dlnbeta_dy_2d = attenuation(ip).dlnbeta_dy_2d;
    subplot(M,N,ip)
    ax = worldmap(lalim, lolim);
    surfacem(xi,yi,beta_tau_2d);
    quiverm(xi,yi,dlnbeta_dy_2d,dlnbeta_dx_2d,'-k')
%     if ~isempty(stlas) 
%         plotm(stlas,stlos,'v');
%     end
    title([num2str(periods(ip)),' s'],'fontsize',15)
    cb = colorbar;
    ylabel(cb,'\nabla\beta/\beta');
%     caxis([-1e-4 4e-4]);
    colormap(flip(seiscmap))
end

% Plot 2D maps of Beta
figure(64); clf; set(gcf,'position',[146           1         726        1024],'color','w');
N=3; M = floor(length(periods)/N)+1;
for ip = 1:length(attenuation)    
    beta_2d = attenuation(ip).beta_2d;
    dlnbeta_dx_2d = attenuation(ip).dlnbeta_dx_2d;
    dlnbeta_dy_2d = attenuation(ip).dlnbeta_dy_2d;
    subplot(M,N,ip)
    ax = worldmap(lalim, lolim);
    surfacem(xi,yi,beta_2d);
    quiverm(xi,yi,dlnbeta_dy_2d,dlnbeta_dx_2d,'-k')
%     if ~isempty(stlas) 
%         plotm(stlas,stlos,'v');
%     end
    title([num2str(periods(ip)),' s'],'fontsize',15)
    cb = colorbar;
    ylabel(cb,'\beta');
    caxis( abs([-1 1]+(max(abs(beta_2d(:)))-1)) );
    colormap(flip(seiscmap))
end

% Re-estimate beta gradient
figure(65); clf; set(gcf,'position',[146           1         726        1024],'color','w');
N=3; M = floor(length(periods)/N)+1;
for ip = 1:length(attenuation)  
    lnbeta = log(attenuation(ip).beta_2d);
    [aspect,slope,dlnbeta_dy_2d,dlnbeta_dx_2d] = gradientm(xi,yi,lnbeta);
    dlnbeta_dy_2d = dlnbeta_dy_2d * 1000;
    dlnbeta_dx_2d = dlnbeta_dx_2d * 1000;
    beta_tau_2d = sqrt(dlnbeta_dx_2d.^2 + dlnbeta_dy_2d.^2);
    subplot(M,N,ip)
    ax = worldmap(lalim, lolim);
%     surfacem(xi,yi,exp(lnbeta));
    surfacem(xi,yi,beta_tau_2d);
    quiverm(xi,yi,dlnbeta_dy_2d,dlnbeta_dx_2d,'-k')
%     if ~isempty(stlas) 
%         plotm(stlas,stlos,'v');
%     end
    title([num2str(periods(ip)),' s'],'fontsize',15)
    cb = colorbar;
    ylabel(cb,'\nabla\beta/\beta Predicted');
%     caxis([-1e-4 4e-4]);
    colormap(flip(seiscmap))
end

% Plot Chi2 of beta fit
figure(66); clf; set(gcf,'color','w');
% subplot(2,1,1);
% for ip = 1:length(periods)
%     plot(periods(ip),receiver(ip).std_pair,'ob','linewidth',2); hold on;
% end
% xlabel('Period (s)');
% ylabel('\sigma log residuals');
% set(gca,'linewidth',1.5,'fontsize',15);
subplot(2,1,2);
plot(periods,[attenuation(:).beta_2d_chi2],'-or','linewidth',2); hold on;
xlabel('Period (s)');
ylabel('\chi^2 misfit');
set(gca,'linewidth',1.5,'fontsize',15);

%% Plot amplitude decay and focusing terms
figure(44); clf;
set(gcf,'Position',[616    71   850   947]);
for ip = 1:length(attenuation)
    subplot(M,N,ip)
    amp_decay_map_evs = squeeze(nanmean(nanmean(attenuation(ip).amp_decay_map,1),2));
    tp_focus_map_evs = squeeze(nanmean(nanmean(attenuation(ip).tp_focus_map,1),2));
    azi_evs = squeeze(nanmean(nanmean(attenuation(ip).azi,1),2));
    plot(attenuation(ip).azi(:),attenuation(ip).amp_decay_map(:),'.b'); hold on;
    plot(attenuation(ip).azi(:),attenuation(ip).tp_focus_map(:),'.r');
%     plot(azi_evs,amp_decay_map_evs,'.b'); hold on;
%     plot(azi_evs,tp_focus_map_evs,'.r');
    title([num2str(attenuation(ip).period),' s'])
    xlabel('azimuth');
    ylabel('amplitude term');
    if ip == 1
        lg = legend({'Amp decay','Focusing'});
        lg.Position(2) = lg.Position(2)+0.07;
    end
    ylim([-2e-4 2e-4]);
end

figure(45); clf;
set(gcf,'Position',[616    71   850   947]);
for ip = 1:length(attenuation)
    subplot(M,N,ip)
    plot(attenuation(ip).tp_focus_map(:),attenuation(ip).amp_decay_map(:),'.b'); hold on;
%     plot(squeeze(azi(9,9,:)),squeeze(amp_term(9,9,:)),'or'); hold on;
%     plot(azi(Ibad),amp_term(Ibad),'o');
%     plot(azi_bin(:),amp_bin(:),'og');
    title([num2str(attenuation(ip).period),' s'])
    ylabel('Amplitude Decay');
    xlabel('Focusing');
end

%% Plot amplitude decay and focusing terms with distance

figure(50); clf;
set(gcf,'Position',[616    71   850   947]);
for ip = 1:length(attenuation)
    subplot(M,N,ip)
    amp_decay_map_evs = squeeze(nanmean(nanmean(attenuation(ip).amp_decay_map,1),2));
    tp_focus_map_evs = squeeze(nanmean(nanmean(attenuation(ip).tp_focus_map,1),2));
    dist_evs = squeeze(nanmean(nanmean(attenuation(ip).dist_map,1),2));
    plot(attenuation(ip).dist_map(:),attenuation(ip).amp_decay_map(:),'.b'); hold on;
    plot(attenuation(ip).dist_map(:),attenuation(ip).tp_focus_map(:),'.r');
%     plot(azi_evs,amp_decay_map_evs,'.b'); hold on;
%     plot(azi_evs,tp_focus_map_evs,'.r');
    
    % Analytical predictions
    x = min(dist_evs):max(dist_evs);
    plot(x,cosd(x)./(sind(x)*6371*nanmean(avgphv(ip).GV_cor(:))),'-m','linewidth',2);
    plot(x,-cosd(x)./(sind(x)*6371*nanmean(avgphv(ip).GV_cor(:))),'-c','linewidth',2);
    plot(x,-cosd(x)./(sind(x)*6371*nanmean(avgphv(ip).GV_cor(:)))-2*alphas(ip)./nanmean(avgphv(ip).GV_cor(:)),'--c','linewidth',2);
    
    title([num2str(attenuation(ip).period),' s'])
    xlabel('distance (deg)');
    ylabel('amplitude term');
    if ip == 1
        lg = legend({'Amp decay','Focusing'});
        lg.Position(2) = lg.Position(2)+0.07;
    end
    ylim([-2e-4 2e-4]);
end

% Sum of terms 
figure(51); clf;
set(gcf,'Position',[616    71   850   947]);
sgtitle('$\frac{c}{2}(\frac{2\nabla{A}\cdot\nabla\tau}{A} + \nabla^2\tau)$','interpreter','latex','fontsize',30);
for ip = 1:length(attenuation)
    subplot(M,N,ip)
    amp_term_map_evs = squeeze(nanmean(nanmean(attenuation(ip).amp_term_2d,1),2));
    dist_evs = squeeze(nanmean(nanmean(attenuation(ip).dist_map,1),2));
    plot(attenuation(ip).dist_map(:),attenuation(ip).amp_term_2d(:),'.b'); hold on;
%     plot(dist_evs,amp_term_map_evs,'.r'); hold on;
%     plot(azi_evs,tp_focus_map_evs,'.r');
    title([num2str(attenuation(ip).period),' s'])
    xlabel('distance (deg)');
    ylabel('amplitude term');
    ylim([-5e-4 2e-4]);
end

%%
figure(61); clf;
set(gcf,'Position',[616    71   850   947]);
sgtitle('$\nabla^2\tau$','interpreter','latex','fontsize',30);
for ip = 1:length(attenuation)
    subplot(M,N,ip)
    tp_focus_map_evs = squeeze(nanmean(nanmean(attenuation(ip).tp_focus_map,1),2));
    dist_evs = squeeze(nanmean(nanmean(attenuation(ip).dist_map,1),2));
%     plot(attenuation(ip).dist_map(:),attenuation(ip).tp_focus_map(:),'.b'); hold on;
    plot(dist_evs,tp_focus_map_evs,'.r'); hold on;
%     plot(azi_evs,tp_focus_map_evs,'.r');
    title([num2str(attenuation(ip).period),' s'])
    xlabel('distance (deg)');
    ylabel('\nabla^2\tau');
    ylim([-0.5e-4 1e-4]);
    
    x = min(dist_evs):max(dist_evs);
    plot(x,cosd(x)./(sind(x)*6371*nanmean(avgphv(ip).GV_cor(:))),'-k','linewidth',2);
% plot(x,-cosd(x)./(sind(x)*6371*nanmean(avgphv(ip).GV_cor(:))),'-c','linewidth',2);
% plot(x,-cosd(x)./(sind(x)*6371*nanmean(avgphv(ip).GV_cor(:)))-2*alphas(ip)./nanmean(avgphv(ip).GV_cor(:)),'--c','linewidth',2);
    set(gca,'fontsize',15)
end

figure(62); clf;
set(gcf,'Position',[616    71   850   947]);
sgtitle('$\frac{2\nabla{A}\cdot\nabla\tau}{A}$','interpreter','latex','fontsize',30);
for ip = 1:length(attenuation)
    subplot(M,N,ip)
    amp_decay_map_evs = squeeze(nanmean(nanmean(attenuation(ip).amp_decay_map,1),2));
    dist_evs = squeeze(nanmean(nanmean(attenuation(ip).dist_map,1),2));
%     plot(attenuation(ip).dist_map(:),attenuation(ip).tp_focus_map(:),'.b'); hold on;
    plot(dist_evs,amp_decay_map_evs,'.r'); hold on;
%     plot(azi_evs,tp_focus_map_evs,'.r');
    title([num2str(attenuation(ip).period),' s'])
    xlabel('distance (deg)');
    ylabel('$\frac{2\nabla{A}\cdot\nabla\tau}{A}$','interpreter','latex');
    ylim([-1.2e-4 0.5e-4]);
    
    x = min(dist_evs):max(dist_evs);
%     plot(x,cosd(x)./(sind(x)*6371*nanmean(avgphv(ip).GV_cor(:))),'-k','linewidth',2);
    plot(x,-cosd(x)./(sind(x)*6371*nanmean(avgphv(ip).GV_cor(:))),'-k','linewidth',2);
    plot(x,-cosd(x)./(sind(x)*6371*nanmean(avgphv(ip).GV_cor(:)))-2*alphas(ip)./nanmean(avgphv(ip).GV_cor(:)),'--k','linewidth',2);
    set(gca,'fontsize',15)
end

%% Auxiliary figures
if is_figures_aux
    %% Plot terms
    figure(46); clf;
    set(gcf,'Position',[616    71   850   947]);
    sgtitle('1 / \nabla\tau','fontweight','bold','fontsize',18);
    for ip = 1:length(attenuation)

        [~,idx] = min(abs(mode.T-periods(ip)));
        phv_tru = mode.phv(idx);

        subplot(M,N,ip)
        r = 0.02;
        tp_grad_map_evs = squeeze(nanmean(nanmean(attenuation(ip).tp_grad_map,1),2));
        azi_evs = squeeze(nanmean(nanmean(attenuation(ip).azi,1),2));
        phvavg = 1./nanmean(attenuation(ip).tp_grad_map(:));
        plot(attenuation(ip).azi(:),1./attenuation(ip).tp_grad_map(:),'.b'); hold on;
        plot(azi_evs,1./tp_grad_map_evs,'.r'); hold on;
        plot([0 360],[phvavg phvavg],'--g','linewidth',2);
    %     plot(azi_evs,tp_grad_map_evs,'.r');
        plot([0 360],[phv_tru phv_tru],'-g','linewidth',2);
        title([num2str(attenuation(ip).period),' s'])
        xlabel('azimuth');
        ylabel('phase velocity');
        if ip == 1
            lg = legend({'\nabla \tau^{-1}'});
            lg.Position(2) = lg.Position(2)+0.07;
        end
        ylim([1-r 1+r]*phvavg);
    end

    figure(47); clf;
    set(gcf,'Position',[616    71   850   947]);
    sgtitle('\nabla A','fontweight','bold','fontsize',18);
    for ip = 1:length(attenuation)
        subplot(M,N,ip)
        amp_grad_map_evs = squeeze(nanmean(nanmean(attenuation(ip).amp_grad_map,1),2));
        azi_evs = squeeze(nanmean(nanmean(attenuation(ip).azi,1),2));
        plot(attenuation(ip).azi(:),(attenuation(ip).amp_grad_map(:)),'.b'); hold on;
        plot(azi_evs,(amp_grad_map_evs),'.r'); hold on;
        title([num2str(attenuation(ip).period),' s'])
        xlabel('azimuth');
        ylabel('amplitude term');
        if ip == 1
            lg = legend({'\nabla A'});
            lg.Position(2) = lg.Position(2)+0.07;
        end
    %     ylim([-2e-4 2e-4]);
    end

    figure(48); clf;
    set(gcf,'Position',[616    71   850   947]);
    sgtitle('\nabla A / A','fontweight','bold','fontsize',18);
    for ip = 1:length(attenuation)
        subplot(M,N,ip)
        amp_grad_norm_map_evs = squeeze(nanmean(nanmean(attenuation(ip).amp_grad_norm_map,1),2));
        azi_evs = squeeze(nanmean(nanmean(attenuation(ip).azi,1),2));
        plot(attenuation(ip).azi(:),(attenuation(ip).amp_grad_norm_map(:)),'.b'); hold on;
        plot(azi_evs,amp_grad_norm_map_evs,'.r'); hold on;
        title([num2str(attenuation(ip).period),' s'])
        xlabel('azimuth');
        ylabel('amplitude term');
        if ip == 1
            lg = legend({'\nabla A / A'});
            lg.Position(2) = lg.Position(2)+0.07;
        end
    %     ylim([-2e-4 2e-4]);
    end

    figure(49); clf;
    set(gcf,'Position',[616    71   850   947]);
    sgtitle('\nabla A \cdot \nabla \tau','fontweight','bold','fontsize',18);
    for ip = 1:length(attenuation)
        subplot(M,N,ip)
        ampgrad_dot_tpgrad = squeeze(nanmean(nanmean(attenuation(ip).ampgrad_dot_tpgrad,1),2));
        azi_evs = squeeze(nanmean(nanmean(attenuation(ip).azi,1),2));
        plot(attenuation(ip).azi(:),(attenuation(ip).ampgrad_dot_tpgrad(:)),'.b'); hold on;
        plot(azi_evs,ampgrad_dot_tpgrad,'.r'); hold on;
        title([num2str(attenuation(ip).period),' s'])
        xlabel('azimuth');
        ylabel('amplitude term');
        if ip == 1
            lg = legend({'\nabla A \cdot \nabla \tau'});
            lg.Position(2) = lg.Position(2)+0.07;
        end
    %     ylim([-2e-4 2e-4]);
    end

    %% Plot amp term with distance

    figure(52); clf;
    set(gcf,'Position',[616    71   850   947]);
    sgtitle('$|\nabla A|$','interpreter','latex','fontsize',30);
    for ip = 1:length(attenuation)
        subplot(M,N,ip)
        amp_grad_map_evs = squeeze(nanmean(nanmean(attenuation(ip).amp_grad_map,1),2));
        dist_evs = squeeze(nanmean(nanmean(attenuation(ip).dist_map,1),2));
    %     plot(attenuation(ip).dist_map(:),attenuation(ip).amp_grad_map(:),'.b'); hold on;
        plot(dist_evs,amp_grad_map_evs,'.r'); hold on;
    %     plot(azi_evs,tp_focus_map_evs,'.r');
        title([num2str(attenuation(ip).period),' s'])
        xlabel('distance (deg)');
        ylabel('|\nabla A|');
    %     ylim([-5e-4 2e-4]);
        ylim([-0.5e-8 0.5e-8]);
    end

    figure(53); clf;
    set(gcf,'Position',[616    71   850   947]);
    sgtitle('$|\nabla \tau|$','interpreter','latex','fontsize',30);
    for ip = 1:length(attenuation)
        subplot(M,N,ip)
        tp_grad_map_evs = squeeze(nanmean(nanmean(attenuation(ip).tp_grad_map,1),2));
        dist_evs = squeeze(nanmean(nanmean(attenuation(ip).dist_map,1),2));
    %     plot(attenuation(ip).dist_map(:),attenuation(ip).tp_grad_map(:),'.b'); hold on;
        plot(dist_evs,tp_grad_map_evs,'.r'); hold on;
    %     plot(azi_evs,tp_focus_map_evs,'.r');
        title([num2str(attenuation(ip).period),' s'])
        xlabel('distance (deg)');
        ylabel('|\nabla \tau|');
    %     ylim([-5e-4 2e-4]);
        ylim(nanmean(tp_grad_map_evs(:))*[0.98 1.02]);
    end

    figure(54); clf;
    set(gcf,'Position',[616    71   850   947]);
    sgtitle('$A$','interpreter','latex','fontsize',30);
    for ip = 1:length(attenuation)
        subplot(M,N,ip)
        amp_map_evs = squeeze(nanmean(nanmean(attenuation(ip).amp_map,1),2));
        dist_evs = squeeze(nanmean(nanmean(attenuation(ip).dist_map,1),2));
    %     plot(attenuation(ip).dist_map(:),attenuation(ip).amp_map(:),'.b'); hold on;
        plot(dist_evs,amp_map_evs,'.r'); hold on;
    %     plot(azi_evs,tp_focus_map_evs,'.r');
        title([num2str(attenuation(ip).period),' s'])
        xlabel('distance (deg)');
        ylabel('A');
    %     ylim([-5e-4 2e-4]);
        set(gca,'yscale','log');
    end

    figure(55); clf;
    set(gcf,'Position',[616    71   850   947]);
    sgtitle('$\nabla A \cdot \nabla \tau$','interpreter','latex','fontsize',30);
    for ip = 1:length(attenuation)
        subplot(M,N,ip)
        ampgrad_dot_tpgrad_evs = squeeze(nanmean(nanmean(attenuation(ip).ampgrad_dot_tpgrad,1),2));
        dist_evs = squeeze(nanmean(nanmean(attenuation(ip).dist_map,1),2));
        plot(attenuation(ip).dist_map(:),attenuation(ip).ampgrad_dot_tpgrad(:),'.b'); hold on;
        plot(dist_evs,ampgrad_dot_tpgrad_evs,'.r'); hold on;
    %     plot(azi_evs,tp_focus_map_evs,'.r');
        title([num2str(attenuation(ip).period),' s'])
        xlabel('distance (deg)');
        ylabel('\nabla A \cdot \nabla \tau');
    %     ylim([-5e-4 2e-4]);
    end

    figure(56); clf;
    set(gcf,'Position',[616    71   850   947]);
    sgtitle('$\nabla A \cdot \nabla \tau$','interpreter','latex','fontsize',30);
    for ip = 1:length(attenuation)
        subplot(M,N,ip)
        ampgrad_dot_tpgrad_ampnorm_evs = squeeze(nanmean(nanmean(attenuation(ip).ampgrad_dot_tpgrad_ampnorm,1),2));
        dist_evs = squeeze(nanmean(nanmean(attenuation(ip).dist_map,1),2));
        plot(attenuation(ip).dist_map(:),attenuation(ip).ampgrad_dot_tpgrad_ampnorm(:),'.b'); hold on;
        plot(dist_evs,ampgrad_dot_tpgrad_ampnorm_evs,'.r'); hold on;
    %     plot(azi_evs,tp_focus_map_evs,'.r');
        title([num2str(attenuation(ip).period),' s'])
        xlabel('distance (deg)');
        ylabel('(\nabla{A} \cdot \nabla\tau) / A');
    %     ylim([-5e-4 2e-4]);
    end

    figure(57); clf;
    set(gcf,'Position',[616    71   850   947]);
    sgtitle('$\frac{|\nabla A|}{A}$','interpreter','latex','fontsize',30);
    for ip = 1:length(attenuation)
        subplot(M,N,ip)
        ampgrad_ampnorm_evs = squeeze(nanmean(nanmean(attenuation(ip).amp_grad_norm_map,1),2));
        dist_evs = squeeze(nanmean(nanmean(attenuation(ip).dist_map,1),2));
    %     plot(attenuation(ip).dist_map(:),attenuation(ip).amp_grad_norm_map(:),'.b'); hold on;
        plot(dist_evs,ampgrad_ampnorm_evs,'.r'); hold on;
    %     plot(azi_evs,tp_focus_map_evs,'.r');
        title([num2str(attenuation(ip).period),' s'])
        xlabel('distance (deg)');
        ylabel('\nabla{A} / A');
        ylim([0 1e-3]);
    end

    % tp_gradR = tp_gradlat.*cosd(azi_prop) + tp_gradlon.*sind(azi_prop);
    % tp_gradT = -tp_gradlat.*sind(azi_prop) + tp_gradlon.*cosd(azi_prop);
    % amp_gradR = amp_gradlat.*cosd(azi_prop) + amp_gradlon.*sind(azi_prop);
    % amp_gradT = -amp_gradlat.*sind(azi_prop) + amp_gradlon.*cosd(azi_prop);
    % Rotate amp_grad coordinate system to be parallel to propagation direction
    figure(58); clf;
    set(gcf,'Position',[616    71   850   947]);
    sgtitle('$\nabla{A_R}$','interpreter','latex','fontsize',30);
    for ip = 1:length(attenuation)
        subplot(M,N,ip)
        amp_gradR_evs = squeeze(nanmean(nanmean(attenuation(ip).amp_gradR_map,1),2));
        dist_evs = squeeze(nanmean(nanmean(attenuation(ip).dist_map,1),2));
    %     plot(attenuation(ip).dist_map(:),attenuation(ip).amp_gradR_map(:),'.b'); hold on;
        plot(dist_evs,amp_gradR_evs,'.r'); hold on;
    %     plot(azi_evs,tp_focus_map_evs,'.r');
        title([num2str(attenuation(ip).period),' s'])
        xlabel('distance (deg)');
        ylabel('\nabla{A_R}');
        ylim([-0.5e-8 0.5e-8]);
    end

    % Rotate amp_grad coordinate system to be parallel to propagation direction
    figure(59); clf;
    set(gcf,'Position',[616    71   850   947]);
    sgtitle('$\frac{\nabla{A_R}}{A}$','interpreter','latex','fontsize',30);
    for ip = 1:length(attenuation)
        subplot(M,N,ip)
    %     amp_gradR_ampnorm = amp_gradR ./ (nanmean(nanmean(attenuation(ip).amp_map,1),2));
        amp_gradR_ampnorm_evs = squeeze(nanmean(nanmean(attenuation(ip).amp_gradR_ampnorm_map,1),2));
        dist_evs = squeeze(nanmean(nanmean(attenuation(ip).dist_map,1),2));
    %     plot(attenuation(ip).dist_map(:),attenuation(ip).amp_gradR_ampnorm_map(:),'.b'); hold on;
        plot(dist_evs,amp_gradR_ampnorm_evs,'.r'); hold on;
    %     plot(azi_evs,tp_focus_map_evs,'.r');
        title([num2str(attenuation(ip).period),' s'])
        xlabel('distance (deg)');
        ylabel('\nabla{A_R} / A');
        ylim([-5e-4 2e-4]);
    end

    % Rotate amp_grad coordinate system to be parallel to propagation direction
    figure(60); clf;
    set(gcf,'Position',[616    71   850   947]);
    sgtitle('$\frac{\nabla{A_R}}{A}|\nabla\tau| \approx \frac{\nabla{A}\cdot\nabla\tau}{A}$','interpreter','latex','fontsize',30);
    for ip = 1:length(attenuation)
        subplot(M,N,ip)
        ampgradR_ampnorm_dot_tpgrad_evs = squeeze(nanmean(nanmean(attenuation(ip).ampgradR_ampnorm_dot_tpgrad,1),2));
        dist_evs = squeeze(nanmean(nanmean(attenuation(ip).dist_map,1),2));
        plot(attenuation(ip).dist_map(:),attenuation(ip).ampgradR_ampnorm_dot_tpgrad(:),'.b'); hold on;
        plot(dist_evs,ampgradR_ampnorm_dot_tpgrad_evs,'.r'); hold on;
    %     plot(azi_evs,tp_focus_map_evs,'.r');
        title([num2str(attenuation(ip).period),' s'])
        xlabel('distance (deg)');
        ylabel('(\nabla{A_R} / A)|\nabla\tau|');
        ylim([-5e-4 2e-4]/4);
    end
end

%% Save
matfilename = fullfile(attenuation_path,['attenuation_',parameters.component,'.mat']);
if is_save_mat
    save(matfilename,'attenuation');
    fprintf('\n');
    disp(['Saved to ',matfilename]);
end
% if is_save_amp_fig
%     figdir = [workingdir,'/figs/attenuation/'];
%     if ~exist(figdir)
%         mkdir(figdir);
%     end
%     save2pdf([figdir,eventphv(1).id,'_attenuation_',parameters.component,'.pdf'],39,100);
%     save2pdf([figdir,eventphv(1).id,'_attenuation_',parameters.component,'_StaAmps.pdf'],40,100);
% end