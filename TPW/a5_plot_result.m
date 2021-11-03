clear;
setup_parameters_tpw;
periods = parameters.periods;
workingdir_tpw = parameters_tpw.workingdir;
gridsize = parameters_tpw.gridsize;
lalim = parameters.lalim;
lolim = parameters.lolim;
xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);

%% Load TPW measurements
r = 0.05;
for ip = 1:length(periods)
    period = periods(ip);
    phvfile = [workingdir_tpw,'/','outvel.',num2str(period,'%03d'),'.txt'];
    azifile = [workingdir_tpw,'/','outazi.',num2str(period,'%03d'),'.txt'];
    stacorfile = [workingdir_tpw,'/','outstacor.',num2str(period,'%03d'),'.txt'];
    alphafile = [workingdir_tpw,'/','outalpha.',num2str(period,'%03d'),'.txt'];
    
    vel(ip) = load_phvfile(phvfile,xnode,ynode);
    ani(ip) = load_azianifile(azifile,xnode,ynode);
    atten(ip) = load_alphafile(alphafile);
    
    tpw.periods(ip) = period;
    tpw.phv_1d(ip) = nanmean(vel(ip).phv(:));
    tpw.phv_1d_std(ip) = nanmean(vel(ip).phv_std(:));
    tpw.A2_1d(ip) = nanmean(ani(ip).A2_kms(:)./tpw.phv_1d(ip));
    tpw.A2_1d_std(ip) = nanmean(ani(ip).A2_std_kms(:)./tpw.phv_1d(ip));
    tpw.phi2_1d(ip) = nanmean(ani(ip).phi2(:));
    tpw.phi2_1d_std(ip) = nanmean(ani(ip).phi2_std(:));
    tpw.alpha_1d(ip) = nanmean(atten(ip).alpha(:));
    tpw.alpha_1d_std(ip) = nanmean(atten(ip).alpha_std(:));

end

%% Load ASWMS measurements
comp = parameters.component;
workingdir = [parameters.workingdir];
eventcs_path = [workingdir,'CSmeasure/'];
eikonal_output_path = [workingdir,'eikonal/'];
eikonal_aniso1D_path = workingdir;
attenuation_path = workingdir; %[workingdir,'attenuation/'];
temp = load([eikonal_aniso1D_path,'/eikonal_stack_aniso1D_',comp,'.mat']);
avgphv_aniso = temp.avgphv_aniso;
temp = load([eikonal_aniso1D_path,'/attenuation_',comp,'.mat']);
attenuation = temp.attenuation;
for ip = 1:length(periods)
    period = periods(ip);
    aswms.periods(ip) = period;
    aswms.phv_1d(ip) = nanmean(avgphv_aniso(ip).isophv(:));
    aswms.phv_1d_std(ip) = nanmean(avgphv_aniso(ip).isophv_std(:));
    aswms.A2_1d(ip) = nanmean(avgphv_aniso(ip).aniso_strength(:));
    aswms.A2_1d_std(ip) = nanmean(avgphv_aniso(ip).aniso_strength_std(:));
    aswms.phi2_1d(ip) = nanmean(avgphv_aniso(ip).aniso_azi(:));
    aswms.phi2_1d_std(ip) = nanmean(avgphv_aniso(ip).aniso_azi_std(:));
    aswms.alpha_1d(ip) = attenuation(ip).alpha_1d;
    aswms.alpha_1d_std(ip) = attenuation(ip).alpha_1d_err;
end

%% Phase velocity maps
figure(31); clf
sgtitle('Phase Velocity','fontweight','bold','fontsize',18)
for ip = 1:length(periods)
    period = periods(ip);
    subplot(4,4,ip)
    ax = worldmap(lalim,lolim);
    set(ax, 'Visible', 'off')
    h1=surfacem(vel(ip).lat,vel(ip).lon,vel(ip).phv);
    % drawpng
    caxis(nanmean(vel(ip).phv(:))*[1-r 1+r]);
    colorbar
    load seiscmap
    colormap(seiscmap);
    title([num2str(period),' s'],'fontsize',16);
end

figure(32); clf
sgtitle('Phase Velocity Uncertainty','fontweight','bold','fontsize',18)
for ip = 1:length(periods)
    period = periods(ip);
    subplot(4,4,ip)
    ax = worldmap(lalim,lolim);
    set(ax, 'Visible', 'off')
    h1=surfacem(vel(ip).lat,vel(ip).lon,vel(ip).phv_std);
    % drawpng
    caxis([0 nanmedian(vel(ip).phv_std(:))*2]);
    colorbar
    load seiscmap
    colormap(seiscmap);
    title([num2str(period),' s'],'fontsize',16);
end

%% Plot 1-D Values
path2qfile = '../qfiles/pa5_5km.s0to66.q';
if exist(path2qfile,'file')==2
    mineos = readMINEOS_qfile(path2qfile,0);
    phv_mineos = interp1(mineos.T,mineos.phv,tpw.periods);
    alpha_mineos = mineos.wrad ./ (2*mineos.grv) ./ mineos.q;
    alpha_mineos = interp1(mineos.T,alpha_mineos,tpw.periods);
end

figure(34); clf;
set(gcf,'position',[616          13         560        1005]);

subplot(4,1,1); hold on;
if exist(path2qfile,'file')==2
    plot(tpw.periods,phv_mineos,'-','color',[0.8 0.8 0.8],'linewidth',2);
end
plot(parameters.periods,parameters_tpw.refphv,'og','linewidth',2);
errorbar(aswms.periods,aswms.phv_1d,aswms.phv_1d_std,'o-r','linewidth',2);
errorbar(tpw.periods,tpw.phv_1d,tpw.phv_1d_std,'o-b','linewidth',2);  
xlabel('Period (s)');
ylabel('Phase Velocity (km/s)');
legend('Mineos','Starting','ASWMS','TPW','location','southeast');
set(gca,'linewidth',1.5,'fontsize',15,'box','on');

subplot(4,1,2); hold on;
plot([min(periods) max(periods)],[0 0],'-','color',[0.8 0.8 0.8],'linewidth',2);
errorbar(aswms.periods,aswms.A2_1d*200,aswms.A2_1d_std*100,'o-r','linewidth',2);
errorbar(tpw.periods,tpw.A2_1d*200,tpw.A2_1d_std*100,'o-b','linewidth',2);
xlabel('Period (s)');
ylabel('Peak-to-peak Anisotropy (%)');
% legend('1-D avg.','True','location','southeast');
set(gca,'linewidth',1.5,'fontsize',15,'box','on');

subplot(4,1,3);  hold on;
errorbar(aswms.periods,aswms.phi2_1d,aswms.phi2_1d_std,'o-r','linewidth',2);
errorbar(tpw.periods,tpw.phi2_1d,tpw.phi2_1d_std,'o-b','linewidth',2);
xlabel('Period (s)');
ylabel('Fast Azimuth (\circ)');
set(gca,'linewidth',1.5,'fontsize',15,'box','on');
% plot([min(periods) max(periods)],78*[1 1],'--k','linewidth',1);
% plot([min(periods) max(periods)],118*[1 1],'--k','linewidth',1);

subplot(4,1,4); hold on;
if exist(path2qfile,'file')==2
    plot(tpw.periods,alpha_mineos,'-','color',[0.8 0.8 0.8],'linewidth',2);
end
errorbar(aswms.periods,aswms.alpha_1d,aswms.alpha_1d_std,'o-r','linewidth',2);
errorbar(tpw.periods,tpw.alpha_1d,tpw.alpha_1d_std,'o-b','linewidth',2);
xlabel('Period (s)');
ylabel('\alpha (km^{-1})');
% legend('1-D avg.','True','location','southeast');
set(gca,'linewidth',1.5,'fontsize',15,'box','on');
