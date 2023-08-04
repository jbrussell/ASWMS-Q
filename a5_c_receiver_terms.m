% Estimate station amplification terms following Eddy & Ekstrom (2014)
% EPSL. 
%
% The script reads single-station amplitude data from the
% auto-correlation measurements. The log amplitudes are then differenced
% between each nearby station pair for each event and frequency and finally
% averaged. The average values are then inverted for individual single
% station receiver terms. As there is no information about absolute
% amplitude in the differential amplitude dataset, we enforce the sum of
% log amplitudes to equal zero.
%
% Corrections to single-station amplitude measurements is then
% Amp_corr = Amp / Amp_rec
% where Amp_rec is the receiver term
%
% github.com/jbrussell
% 2021-05

clear;
setup_parameters

% max_sta_dist = 150; % [km] maximum separation allowed between station pairs
% is_azibin = 1; % bin data by propagation azimuth?
% deg_bins = 15; % [deg] size of azimuthal bins in degrees
% avg_type = 'median'; % 'median'; 'mean'
% min_Mw = 5.5; %6.0;

isfigure = 1;
is_save_mat = 1;

max_sta_dist = parameters.max_sta_dist; % [km] maximum separation allowed between station pairs
is_azibin = parameters.is_azibin; % bin data by propagation azimuth?
deg_bins = parameters.deg_bins; % [deg] size of azimuthal bins in degrees
avg_type = parameters.avg_type; % 'median'; 'mean'
min_Mw = parameters.min_Mw_rec; %6.0;

% input path and files
workingdir = parameters.workingdir;
eventcs_path = [workingdir,'CSmeasure/'];

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
comp = parameters.component;

eventfiles = dir([eventcs_path,'/*_cs_',parameters.component,'.mat']);

load seiscmap

if exist('badampsta.lst','file')
	badstnms = textread('badampsta.lst','%s');
	disp('Found Bad amplitude stations:')
	for ista = 1:length(badstnms)
	disp(badstnms(ista))
    end
else
    badstnms = {};
end


%% Collect single-station amplitude measurements
for ip = 1:length(periods)
    amps_sta = [];
    for ie = 1:length(eventfiles)
        %for ie = 59
        % read in data for this event
        clear eventphv eventcs;
        load(fullfile(eventcs_path,eventfiles(ie).name));
        eventid = eventcs.id;
        disp([num2str(periods(ip)),' : ',eventid]);
        if eventcs.Mw < min_Mw
            continue
        end

		%% Get amplitude measurements
		% reset the arrays
		clear stlas stlos amps
		stlas = eventcs.stlas;
		stlos = eventcs.stlos;
		stnms = eventcs.stnms;
        dists = eventcs.dists;
        evid = eventcs.id;
		if exist('badstnms','var')
			list_badstaids = find(ismember(eventcs.stnms,badstnms));
		else
			list_badstaids = [];
		end
		for ista = 1:length(eventcs.autocor)
			if eventcs.autocor(ista).exitflag(ip)>0 && ~ismember(ista,list_badstaids)
				amp = eventcs.autocor(ista).amp(ip);
			else
				amp = NaN;
            end
            % change from power spectrum to amplitude
            amp = amp.^.5;
            
            azi = azimuth(eventcs.evla,eventcs.evlo,stlas(ista),stlos(ista),referenceEllipsoid('wgs84'));
            
            sta = stnms{ista};
            amplitudes(ip).(['s_',sta]).amp(ie) = amp;
            amplitudes(ip).(['s_',sta]).logamp(ie) = log(amp);
            amplitudes(ip).(['s_',sta]).dist(ie) = dists(ista);
            amplitudes(ip).(['s_',sta]).azi(ie) = azi;
            amplitudes(ip).(['s_',sta]).evid{ie} = evid;
            amplitudes(ip).(['s_',sta]).stla = stlas(ista);
            amplitudes(ip).(['s_',sta]).stlo = stlos(ista);
            amplitudes(ip).(['s_',sta]).stnm = stnms(ista);
            
            % Fill empties with nan
            Iempty = find(amplitudes(ip).(['s_',sta]).amp==0);
            amplitudes(ip).(['s_',sta]).amp(Iempty) = nan;
            amplitudes(ip).(['s_',sta]).logamp(Iempty) = nan;
            if ~isempty(Iempty)
                amplitudes(ip).(['s_',sta]).evid(Iempty) = {'NaN'};
            end
        end            
    end % loop of events
end  % loop of periods

%% Calculate log amplitude ratios, build G matrix, and invert for single station corrections
for ip = 1:length(amplitudes)
    G = [];
    dlogAmp_avg = [];
    std_err = [];
    std_save = [];
    stla = []; stlo = [];
    ipair = 0;
    stas = fields(amplitudes(ip));
    for ista1 = 1:length(stas)
        sta1 = stas{ista1};
		stla(ista1) = amplitudes(ip).(['s_',sta1]).stla;
        stlo(ista1) = amplitudes(ip).(['s_',sta1]).stlo;
        if ismember(sta1,badstnms)
            continue
        end
        logamps1 = amplitudes(ip).(['s_',sta1]).logamp;
        evids1 = amplitudes(ip).(['s_',sta1]).evid;
        for ista2 = 1:length(stas)
            sta2 = stas{ista2};
            if ismember(sta2,badstnms)
                continue
            end
            if strcmp(sta1,sta2)
                continue
            end
            stadist = vdist(amplitudes(ip).(['s_',sta1]).stla,amplitudes(ip).(['s_',sta1]).stlo,...
                            amplitudes(ip).(['s_',sta2]).stla,amplitudes(ip).(['s_',sta2]).stlo)/1000;
            if stadist > max_sta_dist
                continue
            end
            
            logamps2 = amplitudes(ip).(['s_',sta2]).logamp;
            evids2 = amplitudes(ip).(['s_',sta2]).evid;
            
            % Loop through measurements for station pairs i,j and calculate
            % differential log amplitude
            dlogAmp = [];
            azis = [];
            imeas = 0;
            for iev = 1:length(evids1)
                iev2 = find(strcmp(evids1(iev),evids2));
                if ~strcmp(evids1(iev),evids2(iev2))
                    disp('WRONG EVENT')
                    continue
                end
                if isempty(iev2)  || isnan(logamps1(iev)) || isnan(logamps2(iev2)) 
                    continue
                end
                imeas = imeas + 1;
                dlogAmp(imeas) = logamps1(iev)-logamps2(iev2);
                
                azis(imeas) = amplitudes(ip).(['s_',sta1]).azi(iev);
            end
            if isempty(dlogAmp)
                disp([sta1,'-',sta2,': no measurements for this station pair']);
                continue
            end
            
            ipair = ipair + 1;
            
            % Bin data azimuthally
            edges = (0:deg_bins:360);
            [~,~,loc]=histcounts(azis,edges);
            dlogAmp_bin = accumarray(loc(:),dlogAmp(:))./accumarray(loc(:),1);
            azis_bin = 0.5*(edges(1:end-1)+edges(2:end));
            azis_bin = azis_bin(1:length(dlogAmp_bin));
            
            % Unbinned average
            Ne = length(dlogAmp);
            if strcmp(avg_type,'mean')
                dlogAmp_avg_bin = nanmean(dlogAmp_bin);
                dlogAmp_avg_unbinned = nansum(dlogAmp) / Ne;
            elseif strcmp(avg_type,'median')
                dlogAmp_avg_bin = nanmedian(dlogAmp_bin);
                dlogAmp_avg_unbinned = nanmedian(dlogAmp);
            else
                error('avg_type must be ''mean'' or ''median''');
            end
            
            % Data vector for inversion
            if is_azibin==1
                dlogAmp_avg(ipair,:) = dlogAmp_avg_bin;
            else
                dlogAmp_avg(ipair,:) = dlogAmp_avg_unbinned;
            end
            
            if isfigure
                figure(1); clf;
                subplot(2,1,1);
                histogram(azis,25);
                xlim([0 360]);
                title([sta1,'-',sta2,' ',num2str(periods(ip)),' s'])
                ylabel('Number of events')
                set(gca,'fontsize',15)
                subplot(2,1,2);
                plot(azis,dlogAmp,'ob'); hold on;
                plot(azis_bin,dlogAmp_bin,'or');
                plot([0 360],dlogAmp_avg(ipair,:)*[1 1],'-k','linewidth',2);
                plot([0 360],dlogAmp_avg_unbinned*[1 1],'-b');
                plot([0 360],dlogAmp_avg_bin*[1 1],'-r');
                xlim([0 360]);
                xlabel('Azimuth (deg)')
                ylabel('ln(A_i) - ln(A_j)')
                set(gca,'fontsize',15)
            end
            
            % Weighting term
            std_err(ipair,:) = std(dlogAmp) / sqrt(Ne);
            
            % G matrix
            g_row = zeros(1,length(stas));
            g_row(1,ista1) = 1;
            g_row(1,ista2) = -1;
            G(ipair,:) = g_row;
            
            std_save(ipair,:) = std(dlogAmp);
        end
    end
    
    % Remove columns for stations with no pairs
    ista_nopairs = find(~any(G,1));
    stas_bad = stas(ista_nopairs);
    G(:,ista_nopairs) = [];
    stas(ista_nopairs) = [];
    stla(ista_nopairs) = [];
    stlo(ista_nopairs) = [];
    
    % Remove measurements without enough data
    ibad = find(std_err==0);
    dlogAmp_avg(ibad) = [];
    std_err(ibad) = [];
    std_save(ibad) = [];
    G(ibad,:) = [];
    
    % Add final row to ensure all amplitude terms sum to zero
    G(end+1,:) = ones(1,length(stas));
    dlogAmp_avg(end+1,:) = 0;
    std_err(end+1,:) = mean(std_err)*5;
    
    % Invert for receiver amplitude terms
    W = diag(1./std_err).^2;
    F = W.^(0.5)*G;
    f = W.^(0.5)*dlogAmp_avg;
    logAmp_rec = (F'*F)\F'*f;
    Amp_rec = exp(logAmp_rec);
    
    % Estimate chi2 misfit
    dlogAmp_avg_pre = G * logAmp_rec;
    e = (dlogAmp_avg(1:end-1,:) - dlogAmp_avg_pre(1:end-1,:)) ./ std_save;
    chi2 = (e'*e)/length(dlogAmp_avg(1:end-1));
    
%     amplitudes(ip).dlogAmp_avg;
%     amplitudes(ip).std_err;
%     amplitudes(ip).G = G;
	receiver(ip).period = periods(ip);
    receiver(ip).Amp_rec = Amp_rec;
    receiver(ip).std_pair = std_save;
    receiver(ip).chi2 = chi2;
    receiver(ip).stas = stas;
    receiver(ip).stlas = stla;
    receiver(ip).stlos = stlo;
    receiver(ip).stas_bad = stas_bad;
end

%% Fit a smooth surface to the receiver terms
for ip = 1:length(periods)
    stlas = receiver(ip).stlas;
    stlos = receiver(ip).stlos;
    Amp_rec = receiver(ip).Amp_rec;
    [Amp_rec_map,mesh_xi,mesh_yi]=gridfit_jg_geo(stlas,stlos,Amp_rec,xnode,ynode,...
                                   'smooth',2,'regularizer','del4','solver','normal');
    receiver(ip).Amp_rec_map = Amp_rec_map';
    receiver(ip).xi = xi;
    receiver(ip).yi = yi;
end


if is_save_mat==1
    save([workingdir,'receiver_terms_',comp,'.mat'],'receiver');
end


%% Plot results
if isfigure
    figure(49); clf; set(gcf,'color','w');
    subplot(2,1,1);
    for ip = 1:length(periods)
        plot(periods(ip),receiver(ip).std_pair,'ob','linewidth',2); hold on;
    end
    xlabel('Period (s)');
    ylabel('\sigma log residuals');
    set(gca,'linewidth',1.5,'fontsize',15);
    subplot(2,1,2);
    plot(periods,[receiver(:).chi2],'-or','linewidth',2); hold on;
    xlabel('Period (s)');
    ylabel('\chi^2 misfit');
    set(gca,'linewidth',1.5,'fontsize',15);
    
    %%    
    figure(47); clf;
    set(gcf,'Position',[84           3         744        1022],'color','w');
    N=3; M = floor(length(periods)/N)+1;
    sgtitle('Receiver terms','fontweight','bold','fontsize',18);
    for ip = 1:length(periods)
        stlas = receiver(ip).stlas;
        stlos = receiver(ip).stlos;
        Amp_rec = receiver(ip).Amp_rec;

        subplot(M,N,ip)
        ax = worldmap(lalim, lolim);
        scatterm(stlas,stlos,100,Amp_rec,'v','filled','markeredgecolor',[0 0 0]);
        title([num2str(periods(ip)),' s'],'fontsize',15)
        caxis([0.7 1.3]);
        cb = colorbar;
        colormap(seiscmap)
    end

    %%
    figure(48); clf;
    set(gcf,'Position',[84           3         744        1022],'color','w');
    sgtitle('Map of receiver terms','fontweight','bold','fontsize',18);
    for ip = 1:length(periods)
        stlas = receiver(ip).stlas;
        stlos = receiver(ip).stlos;
        Amp_rec_map = receiver(ip).Amp_rec_map;

        subplot(M,N,ip)
        ax = worldmap(lalim, lolim);
        surfacem(xi,yi,Amp_rec_map); hold on;
        plotm(stlas,stlos,'v','markeredgecolor',[0 0 0]);
        title([num2str(periods(ip)),' s'],'fontsize',15)
        caxis([0.7 1.3]);
        cb = colorbar;
        colormap(seiscmap)
    end
    
    %%
    figdir = [workingdir,'/figs/receiver/'];
    if ~exist(figdir)
        mkdir(figdir);
    end
    save2pdf([figdir,'receiver_',parameters.component,'StaAmps.pdf'],47,100);
    save2pdf([figdir,'receiver_',parameters.component,'_StaAmpMap.pdf'],48,100);
end



