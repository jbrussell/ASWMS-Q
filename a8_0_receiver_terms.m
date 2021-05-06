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
% github.com/jbrussell
% 2021-05

clear;
setup_parameters

max_sta_dist = 150; % [km] maximum separation allowed between station pairs
is_azibin = 1; % bin data by propagation azimuth?
deg_bins = 15; % [deg] size of azimuthal bins in degrees

isfigure = 1;
is_save_mat = 1;

% input path and files
workingdir = parameters.workingdir;
eventcs_path = [workingdir,'CSmeasure/'];
eikonal_data_path = [workingdir,'eikonal/'];
eikonal_stack_file = [workingdir,'eikonal_stack_',parameters.component];

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
comp = parameters.component;

eventfiles = dir([eikonal_data_path,'/*_eikonal_',parameters.component,'.mat']);

load seiscmap

if exist('badampsta.lst','file')
	badstnms = textread('badampsta.lst','%s');
	disp('Found Bad amplitude stations:')
	for ista = 1:length(badstnms)
	disp(badstnms(ista))
	end
end


%% Collect single-station amplitude measurements
for ip = 1:length(periods)
    amps_sta = [];
    for ie = 1:length(eventfiles)
        %for ie = 59
        % read in data for this event
        clear eventphv eventcs;
        load(fullfile(eikonal_data_path,eventfiles(ie).name));
        eventid = eventphv(1).id;
        disp(eventid);
        eventcsfile = [eventcs_path,'/',eventid,'_cs_',parameters.component,'.mat'];
        if exist(eventcsfile,'file')
            load(eventcsfile);
        else
            disp(['Cannot find CS file for ',eventid,', Skipped']);
            continue;
        end
        if length(eventphv) ~= length(eventcs.avgphv)
            disp('Inconsist of period number for CS file and eikonal file');
            continue;
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
			if eventcs.autocor(ista).exitflag(ip)>0
				amp = eventcs.autocor(ista).amp(ip);
			else
				amp = NaN;
            end
            % change from power spectrum to amplitude
            amp = amp.^.5;
            
            azi = azimuth(eventcs.evla,eventcs.evlo,stlas(ista),stlos(ista),referenceEllipsoid('wgs84'));

            amplitudes(ip).(stnms{ista}).amp(ie) = amp;
            amplitudes(ip).(stnms{ista}).logamp(ie) = log(amp);
            amplitudes(ip).(stnms{ista}).dist(ie) = dists(ista);
            amplitudes(ip).(stnms{ista}).azi(ie) = azi;
            amplitudes(ip).(stnms{ista}).evid{ie} = evid;
            amplitudes(ip).(stnms{ista}).stla = stlas(ista);
            amplitudes(ip).(stnms{ista}).stlo = stlos(ista);
            amplitudes(ip).(stnms{ista}).stnm = stnms(ista);
            
            % Fill empties with nan
            Iempty = find(amplitudes(ip).(stnms{ista}).amp==0);
            amplitudes(ip).(stnms{ista}).amp(Iempty) = nan;
            amplitudes(ip).(stnms{ista}).logamp(Iempty) = nan;
            if ~isempty(Iempty)
                amplitudes(ip).(stnms{ista}).evid(Iempty) = {'NaN'};
            end
        end            
    end % loop of events
end  % loop of periods

%% Calculate log amplitude ratios, build G matrix, and invert for single station corrections
for ip = 1:length(amplitudes)
    G = [];
    dlogAmp_avg = [];
    std_err = [];
    stla = []; stlo = [];
    ipair = 0;
    stas = fields(amplitudes(ip));
    for ista1 = 1:length(stas)
        sta1 = stas{ista1};
        stla(ista1) = amplitudes(ip).(sta1).stla;
        stlo(ista1) = amplitudes(ip).(sta1).stlo;
        logamps1 = amplitudes(ip).(sta1).logamp;
        evids1 = amplitudes(ip).(sta1).evid;
        for ista2 = 1:length(stas)
            sta2 = stas{ista2};
            if strcmp(sta1,sta2)
                continue
            end
            stadist = vdist(amplitudes(ip).(sta1).stla,amplitudes(ip).(sta1).stlo,...
                            amplitudes(ip).(sta2).stla,amplitudes(ip).(sta2).stlo)/1000;
            if stadist > max_sta_dist
                continue
            end
            
            logamps2 = amplitudes(ip).(sta2).logamp;
            evids2 = amplitudes(ip).(sta2).evid;
            
            % Loop through measurements for station pairs i,j and calculate
            % differential log amplitude
            dlogAmp = [];
            azis = [];
            imeas = 0;
            for iev = 1:length(evids1)
                iev2 = find(ismember(evids1(iev),evids2));
                if isempty(iev2)  || isnan(logamps1(iev)) || isnan(logamps2(iev2)) 
                    continue
                end
                imeas = imeas + 1;
                dlogAmp(imeas) = logamps1(iev)-logamps2(iev2);
                
                azis(imeas) = amplitudes(ip).(sta1).azi(iev);
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
            dlogAmp_avg_bin = nanmean(dlogAmp_bin);
            
            % Data vector for inversion
            Ne = length(dlogAmp);
            if is_azibin==1
                dlogAmp_avg(ipair,:) = dlogAmp_avg_bin;
            else
                dlogAmp_avg(ipair,:) = nansum(dlogAmp) / Ne;
            end
            
            if isfigure
                figure(1); clf;
                subplot(2,1,1);
                histogram(azis,25);
                xlim([0 360]);
                title([sta1,'-',sta2])
                ylabel('Number of events')
                set(gca,'fontsize',15)
                subplot(2,1,2);
                plot(azis,dlogAmp,'ob'); hold on;
                plot([0 360],dlogAmp_avg(ipair,:)*[1 1],'-b');
                plot(azis_bin,dlogAmp_bin,'or');
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
        end
    end
    % Add final row to ensure all amplitude terms sum to zero
    G(ipair+1,:) = ones(1,length(stas));
    dlogAmp_avg(ipair+1,:) = 0;
    std_err(ipair+1,:) = mean(std_err)*5;
    
    % Invert for receiver amplitude terms
    W = diag(1./std_err);
    F = W.^(0.5)*G;
    f = W.^(0.5)*dlogAmp_avg;
    logAmp_rec = (F'*F)\F'*f;
    Amp_rec = exp(logAmp_rec);
    
%     amplitudes(ip).dlogAmp_avg;
%     amplitudes(ip).std_err;
%     amplitudes(ip).G = G;
    receiver(ip).Amp_rec = Amp_rec;
    receiver(ip).stas = stas;
    receiver(ip).stlas = stla;
    receiver(ip).stlos = stlo;
end

%% Fit a smooth surface to the receiver terms
for ip = 1:length(periods)
    stlas = receiver(ip).stlas;
    stlos = receiver(ip).stlos;
    Amp_rec = receiver(ip).Amp_rec;
    [Amp_rec_map,mesh_xi,mesh_yi]=gridfit_jg(stlas,stlos,Amp_rec,xnode,ynode,...
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
    figure(47); clf;
    set(gcf,'Position',[84           3         744        1022],'color','w');
    N=3; M = floor(length(periods)/N)+1;
    sgtitle('Reciver terms','fontweight','bold','fontsize',18);
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



