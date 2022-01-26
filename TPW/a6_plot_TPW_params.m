clear;
setup_parameters_tpw;
periods = parameters.periods;
workingdir_tpw = parameters_tpw.workingdir;
gridsize = parameters_tpw.gridsize;
lalim = parameters.lalim;
lolim = parameters.lolim;
xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);

periods = periods(1:end-1);

r = 0.05;
clear tpw
for ip = 1:length(periods)
    period = periods(ip);
    summfile = [workingdir_tpw,'/','summar.',num2str(round(period),'%03d'),'.inp'];
    allfile = [workingdir_tpw,'/','all.',num2str(round(period),'%03d')];
    tpw(ip) = load_tpw_params(summfile,allfile);
%     tpw(ip).period = period;
end

%% Plot

figure(35); clf;
set(gcf,'color','w','position',[159           1        1197        1002]);
for ip = 1:length(periods)
    subplot(4,3,ip);
%     sgtitle('Phase');
    hold on
    plot(tpw(ip).rmsamp,tpw(ip).startamp2./tpw(ip).startamp1,'.');
    % scatter(tpw(ip).rmsamp,tpw(ip).startamp2./tpw(ip).startamp1,50,tpw(ip).amp1_0,'filled');
    % cb = colorbar;
    % caxis([0 0.3]);
    % colormap(jet);
    % ylabel(cb,'Excitation A_1 / A_0');
    xlabel('RMS amp');
    ylabel('A_2 / A_1');
%     plot(abs(tpw(ip).wvaz2),tpw(ip).startamp2./tpw(ip).startamp1,'.');
%     plot(tpw(ip).startamp1,tpw(ip).startamp2,'.')
    title([num2str(periods(ip)), 's']);
    set(gca,'fontsize',15,'linewidth',1.5,'box','on');
end

figure(36); clf;
set(gcf,'color','w','position',[159           1        1197        1002]);
for ip = 1:length(periods)
    subplot(4,3,ip);
%     sgtitle('Phase');
    hold on
    dwvaz12 = angdiff(tpw(ip).wvaz1*pi/180,tpw(ip).wvaz2*pi/180)*180/pi;
    plot(abs(dwvaz12),tpw(ip).startamp2./tpw(ip).startamp1,'.');
    % scatter(abs(dwvaz12),tpw(ip).startamp2./tpw(ip).startamp1,50,tpw(ip).amp1_0,'filled');
    % cb = colorbar;
    % caxis([0 0.3]);
    % colormap(jet);
    % ylabel(cb,'Excitation A_1 / A_0');
    xlabel('|\phi_1 - \phi_2|');
    ylabel('A_2 / A_1');
    title([num2str(periods(ip)), 's']);
    set(gca,'fontsize',15,'linewidth',1.5,'box','on');
end

%% Wave Amplitude
figure(37); clf;
set(gcf,'color','w','position',[159           1        1197        1002]);
for ip = 1:length(periods)
    subplot(4,3,ip);
%     sgtitle('Phase');
    hold on
    plot(tpw(ip).startamp1,tpw(ip).startamp2,'o');
    xlabel('A_1');
    ylabel('A_2');
    title([num2str(periods(ip)), 's']);
    set(gca,'fontsize',15,'linewidth',1.5,'box','on');
end

%% Wave Angle
figure(38); clf;
set(gcf,'color','w','position',[159           1        1197        1002]);
for ip = 1:length(periods)
    subplot(4,3,ip);
%     sgtitle('Phase');
    hold on
    plot(tpw(ip).wvaz1,tpw(ip).wvaz2,'o');
    xlabel('\phi_1');
    ylabel('\phi_2');
    title([num2str(periods(ip)), 's']);
    set(gca,'fontsize',15,'linewidth',1.5,'box','on');
end
