%%% 5/16/16 -- JBR
%
%
%
% COLUMNS: 1   2  3   4     5     6     7   8   9
%          R,RHO,VPV,VSV,QKAPPA,QSHEAR,VPH,VSH,ETA
%


clear all;
parameter_FRECHET

CARDID = param.CARDID;
frechetpath = param.frechetpath;

issavemat = 1;

YLIMS = [0 500];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_ncard = [param.CARDPATH,param.CARD];
h_nqmod = [param.CARDPATH,param.CARDID,'.qmod'];

card = read_model_card(h_ncard);
R = card.rad;
RHO = card.rho;
VPV = card.vpv;
VSV = card.vsv;
QKAPPA = card.qkap;
QSHEAR = card.qmu;
VPH = card.vph;
VSH = card.vsh;
eta = card.eta;

%% Plot card
fig1 = figure(1); clf; set(gcf, 'Color', 'w');
set(gcf,'position',[ 15 227 1227 478]);

subplot(1,5,1); hold on; box on;
plot(VPV/1000,6371-R/1000,'b','linewidth',2)
plot(VPH/1000,6371-R/1000,'r','linewidth',2)
set(gca,'YDir','reverse');
ylabel('Depth (km)','fontsize',18);
xlabel('V_P (km/s)','fontsize',18);
legend('V_{PV}','V_{PH}','location','southwest');
ylim(YLIMS);
xlim([5 10]);
set(gca,'fontsize',16);

subplot(1,5,2); hold on; box on;
plot(VSV/1000,6371-R/1000,'b','linewidth',2)
plot(VSH/1000,6371-R/1000,'r','linewidth',2)
set(gca,'YDir','reverse');
% ylabel('Depth (km)','fontsize',18);
xlabel('V_S (km/s)','fontsize',18);
legend('V_{SV}','V_{SH}','location','southwest');
ylim(YLIMS);
xlim([3 6]);
set(gca,'fontsize',16);

subplot(1,5,3);
aniso_CARD = (VSH-VSV)./(VSH+VSV)*2*100;
% plot(aniso_CARD,6371-R/1000,'-r','linewidth',2); hold on;
plot(VSH.^2./VSV.^2,6371-R/1000,'-r','linewidth',2); hold on;
plot(VPH.^2./VPV.^2,6371-R/1000,'-b','linewidth',1); hold on;
plot([0 0],YLIMS,'--k','linewidth',2);
set(gca,'YDir','reverse');
% ylabel('Depth (km)','fontsize',18);
% xlabel('Anisotropy \Delta{V}/V (%)','fontsize',18);
xlabel('Anisotropy','fontsize',18);
ylim(YLIMS);
% xlim([-8 8]);
xlim([0.95 1.2]);
legend('\xi','\phi','location','southeast');
set(gca,'fontsize',16);

subplot(1,5,4);
plot(eta,6371-R/1000,'k','linewidth',2)
set(gca,'YDir','reverse');
% ylabel('Depth (km)','fontsize',18);
xlabel('\eta','fontsize',18);
ylim(YLIMS);
% xlim([0.94 1.01]);
xlim([0.9 1.01]);
set(gca,'fontsize',16);

subplot(1,5,5);
plot(RHO,6371-R/1000,'k','linewidth',2)
set(gca,'YDir','reverse');
% ylabel('Depth (km)','fontsize',18);
xlabel('\rho','fontsize',18);
ylim(YLIMS);
% xlim([0.94 1.01]);
% xlim([0.9 1.01]);
set(gca,'fontsize',16);

save2pdf([frechetpath,CARDID,'_CARD'],fig1,600);

if issavemat
    struct.R = R;
    struct.RHO = RHO;
    struct.VPV = VPV;
    struct.VSV = VSV;
    struct.QKAPPA = QKAPPA;
    struct.QSHEAR = QSHEAR;
    struct.VPH = VPH;
    struct.VSH = VSH;
    struct.eta = eta;
    if ~exist('./mat_cards/')
        mkdir('./mat_cards/');
    end
    save(['./mat_cards/',CARDID,'.mat'],'struct')
end