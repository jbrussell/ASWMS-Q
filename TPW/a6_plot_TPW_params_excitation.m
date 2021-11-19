clear;

min_excitation_ratio = 0.6; % minimum source excitation ratio to consider

setup_parameters_tpw;
comp = parameters.component;
periods = parameters.periods;
workingdir_tpw = parameters_tpw.workingdir;
gridsize = parameters_tpw.gridsize;
lalim = parameters.lalim;
lolim = parameters.lolim;
xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
min_sta_num = parameters_tpw.min_sta_num;

% input path
workingdir = [parameters.workingdir];
eikonal_output_path = [workingdir,'eikonal/'];

% periods = periods(1:end-1);

r = 0.05;
clear tpw
for ip = 1:length(periods)
    period = periods(ip);
    summfile = [workingdir_tpw,'/','summar.',num2str(round(period),'%03d'),'.inp'];
    allfile = [workingdir_tpw,'/','all.',num2str(round(period),'%03d')];
    tpw(ip) = load_tpw_params(summfile,allfile);
%     tpw(ip).period = period;
end

% read in bad station list, if existed
if exist('badsta.lst')
    badstnms = textread('badsta.lst','%s');
    disp('Found Bad stations:')
    disp(badstnms)
end

evfiles = dir([path2ASWMS_output,'/CSmeasure/*.mat']);
for iev = 1:length(evfiles)
    evfile = [evfiles(iev).folder,'/',evfiles(iev).name];
    temp = load(evfile);
    eventcs = temp.eventcs;
    temp = load([eikonal_output_path,'/',eventcs.id,'_eikonal_',comp,'.mat']);
    eventphv = temp.eventphv;
    for ip = 1:length(periods)
        if exist('badstnms','var')
            badstaids = find(ismember(stnms,badstnms));
        else
            badstaids = [];
        end
        tp = eventphv(ip).traveltime;
        goodstaind = find(~isnan(tp));
        if ~isempty(badstaids)
            badind = find(ismember(goodstaind,badstaids));
            goodstaind(badind) = [];
        end
        goodstanum = length(goodstaind);
        if goodstanum <= min_sta_num
            continue
        end
        amp1_0 = [];
        ratio_AmpMax = [];
%         for ista = 1:length(eventcs.source)
        for ii = 1:length(goodstaind)
            ista = goodstaind(ii);
            amp1_0(ii) = eventcs.source(ista).excitation(ip).ratio_Amp1_0;
            ratio_AmpMax(ii) = eventcs.source(ista).excitation(ip).ratio_AmpMax(1);
        end
        if nanmean(ratio_AmpMax) < min_excitation_ratio
            tpw(ip).amp1_0(iev) = nan;
            tpw(ip).ratio_AmpMax(iev) = nan;
        else
            tpw(ip).amp1_0(iev) = nanmean(amp1_0);
            tpw(ip).ratio_AmpMax(iev) = nanmean(ratio_AmpMax);
        end
    end
end

for ip = 1:length(periods)
    Inan = isnan(tpw(ip).ratio_AmpMax);
    tpw(ip).amp1_0(Inan) = [];
    tpw(ip).ratio_AmpMax(Inan) = [];
end

%% Plot

figure(35); clf;
set(gcf,'color','w','position',[159           1        1197        1002]);
for ip = 1:length(periods)
    subplot(4,3,ip);
%     sgtitle('Phase');
    hold on
%     plot(tpw(ip).rmsamp,tpw(ip).startamp2./tpw(ip).startamp1,'.');
    scatter(tpw(ip).rmsamp,tpw(ip).startamp2./tpw(ip).startamp1,50,tpw(ip).amp1_0,'filled');
    cb = colorbar;
    caxis([0 1]);
    colormap(jet);
    ylabel(cb,'Excitation A_1 / A_0');
    xlabel('RMS amp');
    ylabel('A_2 / A_1');
%     plot(abs(tpw(ip).wvaz2),tpw(ip).startamp2./tpw(ip).startamp1,'.');
%     plot(tpw(ip).startamp1,tpw(ip).startamp2,'.')
    title([num2str(periods(ip)), 's']);
    set(gca,'fontsize',15,'linewidth',1.5,'box','on');
end

%%
figure(36); clf;
set(gcf,'color','w','position',[159           1        1197        1002]);
for ip = 1:length(periods)
    subplot(4,3,ip);
%     sgtitle('Phase');
    hold on
    dwvaz12 = angdiff(tpw(ip).wvaz1*pi/180,tpw(ip).wvaz2*pi/180)*180/pi;
%     plot(abs(dwvaz12),tpw(ip).startamp2./tpw(ip).startamp1,'.');
    scatter(abs(dwvaz12),tpw(ip).startamp2./tpw(ip).startamp1,50,tpw(ip).amp1_0,'filled');
    cb = colorbar;
    caxis([0 1]);
    colormap(jet);
    ylabel(cb,'Excitation A_1 / A_0');
    xlabel('|\phi_1 - \phi_2|');
    ylabel('A_2 / A_1');
    title([num2str(periods(ip)), 's']);
    set(gca,'fontsize',15,'linewidth',1.5,'box','on');
end

%%
figure(38); clf;
set(gcf,'color','w','position',[159           1        1197        1002]);
for ip = 1:length(periods)
    subplot(4,3,ip);
%     sgtitle('Phase');
    hold on
    dwvaz12 = angdiff(tpw(ip).wvaz1*pi/180,tpw(ip).wvaz2*pi/180)*180/pi;
%     plot(abs(dwvaz12),tpw(ip).startamp2./tpw(ip).startamp1,'.');
    scatter(tpw(ip).rmsamp,tpw(ip).amp1_0,50,abs(dwvaz12),'filled');
    cb = colorbar;
%     caxis([0 0.3]);
    colormap(jet);
    ylabel(cb,'|\phi_1 - \phi_2|');
    xlabel('RMS amp');
    ylabel('Excitation A_1 / A_0');
    title([num2str(periods(ip)), 's']);
    set(gca,'fontsize',15,'linewidth',1.5,'box','on');
end

%%
figure(39); clf;
set(gcf,'color','w','position',[159           1        1197        1002]);
for ip = 1:length(periods)
    subplot(4,3,ip);
%     sgtitle('Phase');
    hold on
%     plot(tpw(ip).rmsamp,tpw(ip).startamp2./tpw(ip).startamp1,'.');
    scatter(tpw(ip).rmsamp,tpw(ip).startamp2./tpw(ip).startamp1,50,tpw(ip).ratio_AmpMax,'filled');
    cb = colorbar;
    caxis([0 1]);
    colormap(jet);
    ylabel(cb,'A_0/A_{max}');
    xlabel('RMS amp');
    ylabel('A_2 / A_1');
%     plot(abs(tpw(ip).wvaz2),tpw(ip).startamp2./tpw(ip).startamp1,'.');
%     plot(tpw(ip).startamp1,tpw(ip).startamp2,'.')
    title([num2str(periods(ip)), 's']);
    set(gca,'fontsize',15,'linewidth',1.5,'box','on');
end

%%
figure(40); clf;
set(gcf,'color','w','position',[159           1        1197        1002]);
for ip = 1:length(periods)
    subplot(4,3,ip);
%     sgtitle('Phase');
    hold on
    dwvaz12 = angdiff(tpw(ip).wvaz1*pi/180,tpw(ip).wvaz2*pi/180)*180/pi;
%     plot(abs(dwvaz12),tpw(ip).startamp2./tpw(ip).startamp1,'.');
    scatter(abs(dwvaz12),tpw(ip).startamp2./tpw(ip).startamp1,50,tpw(ip).ratio_AmpMax,'filled');
    cb = colorbar;
    caxis([0 1]);
    colormap(jet);
    ylabel(cb,'A_0/A_{max}');
    xlabel('|\phi_1 - \phi_2|');
    ylabel('A_2 / A_1');
    title([num2str(periods(ip)), 's']);
    set(gca,'fontsize',15,'linewidth',1.5,'box','on');
end