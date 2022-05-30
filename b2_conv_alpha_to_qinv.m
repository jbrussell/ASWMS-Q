% Convert estimates of alpha to Q^-1
%
% Q^-1 = 2 * alpha * grv / omega
%
% github.com/jbrussell
% 2022-05

clear;
% setup parameters
setup_parameters

is_save_mat = 1;

is_eikonal_ampgrad_norm = parameters.is_eikonal_ampgrad_norm;

r = 0.05;

workingdir = parameters.workingdir;
eikonal_grv_stack_file = [workingdir,'eikonal_grv_stack_',parameters.component];
attenuation_path = workingdir; %[workingdir,'attenuation/'];
attenuation_file = [attenuation_path,'attenuation_',parameters.component];

% load stacked phase velocity map
temp = load(eikonal_grv_stack_file);
avggrv = temp.avggrv;
% Load attenuation
temp = load(attenuation_file);
attenuation = temp.attenuation;

% set up useful variables
lalim = parameters.lalim;
lolim = parameters.lolim;
gridsize = parameters.gridsize;
xnode = lalim(1):gridsize:lalim(2);
ynode = lolim(1):gridsize:lolim(2);
[xi yi] = ndgrid(xnode,ynode);
Nx=length(xnode); Ny=length(ynode);
periods = parameters.periods;

load seiscmap

%% Do conversions

for ip = 1:length(attenuation)
    
    % Unbinned 1-D sinusoidal fit
    attenuation(ip).qinv_1d = attenuation(ip).alpha_1d*2.*nanmean(avggrv(ip).GV(:)) ./ (2*pi./attenuation(ip).period);
    attenuation(ip).qinv_1d_err = 2./(2*pi./attenuation(ip).period) .* sqrt((attenuation(ip).alpha_1d_err.*nanmean(avggrv(ip).GV(:))).^2 + (attenuation(ip).alpha_1d.*nanmean(avggrv(ip).GV_std(:))).^2);
    
    % Binned 1-D sinusoidal fit
    attenuation(ip).qinv_1d_bin = attenuation(ip).alpha_1d_bin*2.*nanmean(avggrv(ip).GV(:)) ./ (2*pi./attenuation(ip).period);
    attenuation(ip).qinv_1d_bin_err = 2./(2*pi./attenuation(ip).period) .* sqrt((attenuation(ip).alpha_1d_bin_err.*nanmean(avggrv(ip).GV(:))).^2 + (attenuation(ip).alpha_1d_bin.*nanmean(avggrv(ip).GV_std(:))).^2);
    
    % 2-D sinusoidal fit
    attenuation(ip).qinv_2d = attenuation(ip).alpha_2d*2.*avggrv(ip).GV ./ (2*pi./attenuation(ip).period);
    attenuation(ip).qinv_2d_err = 2./(2*pi./attenuation(ip).period) .* sqrt((attenuation(ip).alpha_2d_err.*avggrv(ip).GV).^2 + (attenuation(ip).alpha_2d.*avggrv(ip).GV_std).^2);
    
    % Get average values from 2-D maps
    attenuation(ip).qinv_2d_mean = nanmean(attenuation(ip).qinv_2d(:));
    attenuation(ip).qinv_2d_mean_err = max([nanstd(attenuation(ip).qinv_2d(:)), nanmean(attenuation(ip).qinv_2d_err(:))]);
    % Get alpha from center of array
    latc = mean(xi(~isnan(attenuation(ip).qinv_2d))); % center latitude
    lonc = mean(yi(~isnan(attenuation(ip).qinv_2d))); % center longitude
    [~,ilat] = min(abs(xnode-latc));
    [~,ilon] = min(abs(ynode-lonc));
    qinv_2d_block = attenuation(ip).qinv_2d(ilat+[-1:1],ilon+[-1:1]);
    qinv_2d_block_std = attenuation(ip).qinv_2d_err(ilat+[-1:1],ilon+[-1:1]);
    attenuation(ip).qinv_2d_center = nanmean(qinv_2d_block(:));
    attenuation(ip).qinv_2d_center_err = max([nanstd(qinv_2d_block(:)), nanmean(qinv_2d_block_std(:))]);
    attenuation(ip).latc = latc;
    attenuation(ip).lonc = lonc;
    
    attenuation(ip).grv = avggrv(ip).GV;
    attenuation(ip).grv_std = avggrv(ip).GV_std;

end

%% Save
matfilename = fullfile(attenuation_path,['attenuation_',parameters.component,'.mat']);
if is_save_mat
    save(matfilename,'attenuation');
    fprintf('\n');
    disp(['Saved to ',matfilename]);
end

%% Plot 1D average qinv

% % mode = readMINEOS_qfile('./qfiles/pa5_5km.s0to66.q',0);
% % mode = readMINEOS_qfile('./qfiles/S362ANI_NoMelt.s0to100.q',0);
% mode = readMINEOS_qfile('./qfiles/S362ANI_JdF.s0to66.q',0);
% alpha_MINEOS = mode.wrad ./ (2*mode.grv) ./ mode.q;

% temp = load(['../write_profile/CARDS_MINEOS_s362ani_crust2.0.mat']);
% mat = temp.mat; clear temp;
% mode.T = [mat(:).period];
% region = [];
% for ip = 1:length(mat)
%     omega = 2*pi ./ mode.T(ip);
%     mat(ip).alpha = omega./(2.*mat(ip).grv .* mat(ip).q);
% %     mat(ip).alpha = omega./(2.*4 .* mat(ip).q);
%     region(ip).q2d = interp2(mat(1).lon,mat(1).lat,mat(ip).q,yi,xi);
%     mode.q(ip) = nanmean(region(ip).q2d(:));
% end
% qinv_MINEOS =  1 ./ mode.q;

figure(42); clf; set(gcf,'color','w');
alpha_zhitu = [4.1 7.3 8.2 8.9 6.9]*1e-5;
f_mhz_zhitu = [10 15 20 25 30];
qinv_zhitu = 2.*alpha_zhitu.*4 ./ (2*pi*f_mhz_zhitu/1000);
qinvs = [attenuation(:).qinv_1d];
qinvs_err = [attenuation(:).qinv_1d_err];
qinvs_bin = [attenuation(:).qinv_1d_bin];
qinvs_bin_err = [attenuation(:).qinv_1d_bin_err];
% qinvs_avg = [attenuation(:).qinv_1d_avg];
% qinvs_avg_err = [attenuation(:).qinv_1d_avg_err];
qinvs_2d = [attenuation(:).qinv_2d_mean];
qinvs_2d_err = [attenuation(:).qinv_2d_mean_err];
qinvs_2d_center = [attenuation(:).qinv_2d_center];
qinvs_2d_center_err = [attenuation(:).qinv_2d_center_err];
% plot(mode.T,qinv_MINEOS,'-','color',[0.7 0.7 0.7],'linewidth',5); hold on;
errorbar(periods,qinvs_bin,qinvs_bin_err,'-om'); hold on;
errorbar(periods,qinvs,qinvs_err,'-ok');
% plot(periods,qinvs_avg,'-oc');
% errorbar(periods,qinvs_2d,qinvs_2d_err,'-ob');
plot(periods,qinvs_2d,'-ob');
errorbar(periods,qinvs_2d_center,qinvs_2d_center_err,'-og');
plot(1./f_mhz_zhitu*1000,qinv_zhitu,'xr','linewidth',3,'MarkerSize',8); hold on;
plot(1./f_mhz_zhitu*1000,qinv_zhitu*2,'x','color',[0 0.85 0],'linewidth',3,'MarkerSize',8); hold on;
legend({'True','1D fit (bin)','1D fit','1D mean','2D fit center','Zhitu','Zhitu x 2'},'location','northeastoutside','fontsize',15)
set(gca,'fontsize',15,'linewidth',1.5);
xlabel('Period (s)');
ylabel('Q^{-1}');
xlim([min(periods)-10 max(periods)+10]);

%% Plot 2D maps of qinv
figure(43); clf; set(gcf,'position',[146           1         726        1024],'color','w');
N=3; M = floor(length(periods)/N)+1;
for ip = 1:length(attenuation)    
    qinv_2d = attenuation(ip).qinv_2d;    
    subplot(M,N,ip)
    ax = worldmap(lalim, lolim);
    surfacem(xi,yi,qinv_2d);
%     if ~isempty(stlas) 
%         plotm(stlas,stlos,'v');
%     end
    title([num2str(periods(ip)),' s'],'fontsize',15)
    cb = colorbar;
    ylabel(cb,'Q^{-1}');
    caxis([0 2e-2]);
    colormap(flip(seiscmap))
end

