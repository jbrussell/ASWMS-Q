clear; close;

setup_parameters
workingdir = parameters.workingdir;

path2file = {
             [workingdir,'/eikonal_ani_BHZ.mat'];
             };

% path2file = {
%              './matgsdf_30stas_12_27_SNR10/ORCA_M5.5_detrend_Zcorr_100km_snr3_600km/eikonal_ani_BHZ_200km_cohere0.8_12_27s.mat';
%              './matgsdf_30stas_20_100/ORCA_M5.5_detrend_Zcorr_100km_snr3_600km/eikonal_ani_BHZ_200km_cohere0.8_20_100s.mat';
%              './matgsdf_30stas_70_150/ORCA_M5.5_detrend_Zcorr_100km_snr3_600km/eikonal_ani_BHZ_200km_cohere0.8_70_150s.mat';
%              };
figdir = ['figs_paper/'];

Iper_plot = 1;
lalim = parameters.lalim;
lolim = parameters.lolim;

%%
for ii = 1:length(path2file)
    data(ii) = load(path2file{ii});
end

eventphv_ani = [data(:).eventphv_ani];
temp.fit_azi = [data(:).fit_azi];
temp.fit_azi_bin = [data(:).fit_azi_bin];
temp.fit_azi_bin_res = [data(:).fit_azi_bin_res];

flds = fields(temp.fit_azi(1));
for ii = 1:length(flds)
    fit_azi.(flds{ii}) = [temp.fit_azi(:).(flds{ii})];
end
flds = fields(temp.fit_azi_bin(1));
for ii = 1:length(flds)
    fit_azi_bin.(flds{ii}) = [temp.fit_azi_bin(:).(flds{ii})];
end
flds = fields(temp.fit_azi_bin_res(1));
for ii = 1:length(flds)
    fit_azi_bin_res.(flds{ii}) = [temp.fit_azi_bin_res(:).(flds{ii})];
end

% './matgsdf_30stas_12_27_SNR10/ORCA_M5.5_detrend_Zcorr_100km_snr3_600km/eikonal_ani_BHZ_200km_cohere0.8_12_27s.mat';
% './matgsdf_30stas_20_100/ORCA_M5.5_detrend_Zcorr_100km_snr3_600km/eikonal_ani_BHZ_200km_cohere0.8_20_100s.mat';
% './matgsdf_30stas_70_150/ORCA_M5.5_detrend_Zcorr_100km_snr3_600km/eikonal_ani_BHZ_200km_cohere0.8_70_150s.mat';
% './matgsdf_30stas_DPG_20_50/ORCA_M5.5_detrend_100km_snr3_600km/eikonal_ani_BDH_200km_cohere0.8_20_50s.mat';

%% Plot map
period = eventphv_ani(Iper_plot).period;
isgood = eventphv_ani(Iper_plot).isgood;
evlas = eventphv_ani(Iper_plot).evla(isgood);
evlos = eventphv_ani(Iper_plot).evlo(isgood);
Mag = eventphv_ani(Iper_plot).Mw(isgood);
[~,Ievs] = unique(evlas .* evlos);
Mag = Mag(Ievs);
lats = evlas(Ievs);
lons = evlos(Ievs);

figure(3); clf;
landareas = shaperead('landareas.shp','UseGeoCoords',true);
ax = axesm('eqdazim', 'Frame', 'on', 'Grid', 'off');
box off;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
% setm(ax,'Origin',[mean(lalim),mean(lolim)])
setm(ax,'Origin',[mean(lalim),mean(lolim)],'FLatLimit',[-125 125]+mean(lalim),'FLonLimit',[],'MapLonLimit',[-125 125]+mean(lolim))
geoshow(ax, landareas,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none'); hold on;
plotm(lats,lons,'o','color',[0.4 0 0],'MarkerFaceColor',[0.85 0 0],'MarkerSize',10,'linewidth',1);
plotm(mean(lalim),mean(lolim),'p','color',[0 0.2 0.4],'MarkerFaceColor',[0 0.5 1],'MarkerSize',24,'linewidth',1);
for ii = [30 60 90 120]
    [latc,longc] = scircle1(mean(lalim),mean(lolim),ii);
    plotm(latc,longc,'-k','linewidth',1)
end

% save2pdf([figdir,'events.pdf'],3,500);
