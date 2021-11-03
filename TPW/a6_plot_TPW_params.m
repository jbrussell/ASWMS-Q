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
    summfile = [workingdir_tpw,'/','summar.',num2str(period,'%03d'),'.inp'];
    allfile = [workingdir_tpw,'/','all.',num2str(period,'%03d')];
    tpw(ip) = load_tpw_params(summfile,allfile);
%     tpw(ip).period = period;
end

%% Plot

figure(35); clf;
for ip = 1:length(periods)
    subplot(4,4,ip);
%     sgtitle('Phase');
    hold on
    plot(tpw(ip).rmsamp,tpw(ip).startamp2./tpw(ip).startamp1,'.');
    xlabel('RMS amp');
    ylabel('A_2 / A_1');
%     plot(abs(tpw(ip).wvaz2),tpw(ip).startamp2./tpw(ip).startamp1,'.');
%     plot(tpw(ip).startamp1,tpw(ip).startamp2,'.')
end