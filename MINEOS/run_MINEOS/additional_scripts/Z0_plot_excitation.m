% Plot source excitations
%

clear; close all;

COMP = 'T'; %'Z' 'R' 'T' % Component
MODE = 0;

parameter_FRECHET;
STAPATH = param.STAPATH;
SYNTH_OUT = param.SYNTH_OUT;
EXCITE_OUT = [SYNTH_OUT,'excitation/'];
EVTPATH = param.EVTPATH;

% Load stations
sta = load_stations(STAPATH);

% Load event info
moment = load_moment_tensor(EVTPATH);

%%
figure(1); clf;
for ista = 1:length(sta.staname)
    path2asc = [EXCITE_OUT,'/',sta.staname{ista},'.',COMP,'.excite.asc'];
    [excite] = load_excitation_asc(path2asc,MODE);
    
    subplot(2,2,1);
    plot(excite.per,excite.amp,'linewidth',1); hold on;
    xlabel('Period');
    ylabel('ak');
    xlim([0 100]);
    
    subplot(2,2,2);
    plot(excite.per,excite.ae,'linewidth',1); hold on;
    xlabel('Period');
    ylabel('ae');
    xlim([0 100]);
    
    subplot(2,2,3);
    plot(excite.per,excite.phase,'linewidth',1); hold on;
    xlabel('Period');
    ylabel('phi2');
    xlim([0 100]);
    
    subplot(2,2,4);
    plot(excite.per,excite.phip,'linewidth',1); hold on;
    xlabel('Period');
    ylabel('phip');
    xlim([0 100]);
end

%%
figure(2); clf;
set(gcf,'Position',[20   605   940*1.7   420],'color','w');

branches = [0 1];
period = 84; % [s]
clrs = jet(length(branches));
for ibr = 1:length(branches)
    clear X az A phi
    for ista = 1:length(sta.staname)
        [X(ista),az(ista)] = distance(moment.lat,moment.lon,sta.lats(ista),sta.lons(ista));

        path2asc = [EXCITE_OUT,'/',sta.staname{ista},'.',COMP,'.excite.asc'];
        [excite] = load_excitation_asc(path2asc,branches(ibr));
        
        A(ista) = interp1(excite.per,excite.amp,period);
%         A(ista) = interp1(excite.per,excite.ae,period);
        phi(ista) = interp1(excite.per,excite.phase,period);
        phi(ista) = angdiff(phi(ista),0);
    end
        
    % Plot amplitude radiation pattern
    if ibr == 1
        ax1 = subplot(1,3,1);
        polarAxesHandle = polaraxes('Units',ax1.Units,'Position',ax1.Position);
        delete(ax1);
        ax1 = polarAxesHandle;
        h1(ibr) = polarplot(ax1,az*pi/180,A,'o','color',clrs(ibr,:),'linewidth',1.5); hold on;
        ph = h1(ibr).Parent; %axis handle 
        hold(ph, 'on')
        ax1.Position = [ax1.Position(1) ax1.Position(2)-0.2 ax1.Position(3) ax1.Position(4)*1.5];
    else
        h1(ibr) = polarplot(ph,az*pi/180,A,'o','color',clrs(ibr,:),'linewidth',1.5); hold on;
    end
    legendstr{ibr} = ['n = ' num2str(ibr-1)];
    hlegend=legend(h1,legendstr);
    title(['\textbf{Amp ' num2str(period) 's}'],'interpreter','latex')
    grid on; box on;
    set(ax1,'fontweight','bold')
    set(ax1,'RTickLabel',[],'ThetaZeroLocation','top','ThetaDir','clockwise');
    set(ax1,'fontsize',20)

    % Plot phase variation with azimuth
    if ibr == 1
        ax2 = subplot(1,3,2); hold on;
    end
    h2(ibr) = plot(ax2,az,phi,'o','color',clrs(ibr,:),'linewidth',1.5); hold on;
    plot(ax2,az,phi-2*pi,'o','color',clrs(ibr,:),'linewidth',1.5); hold on;
    plot(ax2,az,phi+2*pi,'o','color',clrs(ibr,:),'linewidth',1.5); hold on;
    hlegend=legend(h2,legendstr);
    ylim(ax2,[-3*pi/2 3*pi/2]);
    xlim(ax2,[0 360]);
    xlabel(ax2,'Azimuth','interpreter','latex')
    ylabel(ax2,'Phase (rad)','interpreter','latex')
    title(['\textbf{Phase ' num2str(period) 's}'],'interpreter','latex')
    set(ax2,'fontsize',20)
    grid on; box on;
    set(ax2,'fontweight','bold')
    
    amps{ibr} = A;
end

subplot(1,3,3);
A1_A0 = amps{2} ./ amps{1};
plot(az,A1_A0,'o','color',clrs(1,:),'linewidth',1.5);
xlim([0 360]);
xlabel('Azimuth','interpreter','latex')
ylabel('Excitation ratio','interpreter','latex')
title(['\textbf{' num2str(period) 's}'],'interpreter','latex')
set(gca,'fontsize',20)
grid on; box on;
set(gca,'fontweight','bold')
legend('n1/n0');

% Focal mecanism
axes('Position',[-0.06 .65 .25 .25])
axis square; axis off;
focalmech([moment.m_rr,moment.m_tt,moment.m_pp,moment.m_rt,moment.m_rp,moment.m_tp],0,0,10,'text',{moment.evid; [num2str(moment.depth_km),' km']})