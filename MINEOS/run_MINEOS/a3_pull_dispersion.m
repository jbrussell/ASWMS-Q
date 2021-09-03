% 10/8/15 -- Josh Russell
%
% Program to plot l-w dispersion, Group velocity dispersion, and Phase
% velocity dispersion from the *.q ascii file.
%
% ~~ RUN MINEOS FIRST ~~
%
% Columns of dat{i}: 
%               nn,ll,w,qq,phi,cv,gv,cvq,Tq,T
%
clear all;

%%%% PARAMETERS %%%%
bw = 0;        % print in black and white? (1=yes, 0=no)
paus = 0;      % pause figure? (1=yes, 0=no)
Tlim = [5 100]; %[2 100]; %[20 150]; %[0 50]; % sec
gvlim = [1 5]; %[1 5]; %[3 5]; %[1 8]; % km/s
cvlim = [3 5]; %[1 5]; %[3 5]; %[1 8]; % km/s
%%%%%%%%%%%%%%%%%%%%%


fig_lw = figure(1); subplot(2,1,2); hold on; box on; set(gcf, 'Color', 'w'); clf;
fig_gv = figure(2); subplot(2,1,2); hold on; box on; set(gcf, 'Color', 'w'); clf;
fig_cv = figure(3); subplot(2,1,2); hold on; box on; set(gcf, 'Color', 'w'); clf;
fig_cvgv = figure(5); hold on; box on; set(gcf, 'Color', 'w'); clf;

parameter_FRECHET;
TABLEPATH = param.TABLEPATH;
CARDID = param.CARDID;
TYPE = param.TYPE;

% N_modes = 4

if strcmp(TYPE,'S') == 1
    TYPEID = param.STYPEID; 
elseif strcmp(TYPE,'T') == 1
    TYPEID = param.TTYPEID;
end

QIN = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.q'];

if bw == 1
    clr = copper(N_modes);
elseif bw == 0
    clr = jet(N_modes);
end
% clr = [0 0 1; 0 0.7 0];
% clr = [0 0 0; 1 0 0];
clr = lines(N_modes);
dat = {};
for i = 1:N_modes
%     com = ['awk ''{ if ($1 ==',num2str(i-1),' && $10 != "") print $0}'' ',QIN];
%     [log3, dat{i}] = system(com);
%     dat{i} = str2num(dat{i});
%     nn =  dat{i}(:,1);
%     ll =  dat{i}(:,2);
%     w =   dat{i}(:,3)/(2*pi)*1000; %convert rad/s ---> mhz
%     qq =  dat{i}(:,4);
%     phi = dat{i}(:,5);
%     cv =  dat{i}(:,6);
%     gv =  dat{i}(:,7);
%     cvq = dat{i}(:,8);
%     Tq =  dat{i}(:,9);
%     T =   dat{i}(:,10);
    
    [mode] = readMINEOS_qfile(i-1);
    ll = mode.l;
    cv = mode.phv;
    w = mode.w;
    cvq = mode.phvq;
    gv = mode.grv;
    Tq = mode.Tq;
    T = mode.T;
    
    % PLOT L-W
    figure(fig_lw); hold on; box on;
    plot(ll(:),w(:),'-o','color',clr(i,:),'linewidth',2,'markersize',0.5);
    t1 = title(sprintf('%s (maxF = %.1f maxL = %d)',CARDID,maxF,maxL),'fontsize',12);
    set(gca,'fontsize',16,'linewidth',2);
    xlabel('Angular Degree l','fontsize',16);
    ylabel('Frequency \omega (mhz)','fontsize',16);
    %xlim([0 175]);
    %ylim([0 20]);
%     ylim([0 100]);
    
    % PLOT gv Dispersion (T-gv)
    figure(fig_gv);hold on; box on;
    plot(T(:),gv(:),'-o','color',clr(i,:),'linewidth',2,'markersize',0.5);
    t2 = title(sprintf('%s (maxF = %.1f maxL = %d)',CARDID,maxF,maxL),'fontsize',12);
    xlabel('Period T (sec)','fontsize',12);
    ylabel('Group Velocity (km/s)','fontsize',12);
    xlim(Tlim);
    ylim(gvlim);
    
    figure(fig_cvgv); subplot(2,1,2); hold on; box on;
    plot(T(:),gv(:),'-o','color',clr(i,:),'linewidth',3,'markersize',0.5);
    t2 = title(sprintf('%s (maxF = %.1f maxL = %d)',CARDID,maxF,maxL),'fontsize',15);
    xlabel('Period T (sec)','fontsize',15);
    ylabel('Group Velocity (km/s)','fontsize',15);
    xlim(Tlim);
    ylim(gvlim);
    set(gca,'fontsize',15,'linewidth',2,'xminortick','on','yminortick','on')
    
    % PLOT q corrected cv dispersion (Tq-cvq)
    figure(fig_cv);hold on; box on;
    plot(Tq(:),cvq(:),'-o','color',clr(i,:),'linewidth',2,'markersize',0.5);
%     t3 = title(sprintf('%s (maxF = %.1f maxL = %d)',CARDID,maxF,maxL),'fontsize',12);
    xlabel('Period T (sec)','fontsize',12);
    ylabel('Phase Velocity (km/s)','fontsize',12);
    xlim(Tlim);
    ylim(cvlim);
    
    figure(fig_cvgv); subplot(2,1,1); hold on; box on;
    plot(Tq(:),cvq(:),'-o','color',clr(i,:),'linewidth',3,'markersize',0.5);
%     t3 = title(sprintf('%s (maxF = %.1f maxL = %d)',CARDID,maxF,maxL),'fontsize',15);
    xlabel('Period T (sec)','fontsize',15);
    ylabel('Phase Velocity (km/s)','fontsize',15);
    xlim(Tlim);
    ylim(cvlim);
    set(gca,'fontsize',15,'linewidth',2,'xminortick','on','yminortick','on')
    
    if paus
        pause;
    end
    
end

lgd = {};
for i = 1:N_modes
    lgd{i} = [TYPE,num2str(i-1)];
end
figure(fig_lw);
legend(lgd,'location','southeastoutside');
figure(fig_gv);
legend(lgd,'location','northeastoutside');
figure(fig_cv);
legend(lgd,'location','northeastoutside');
figure(fig_cvgv);
legend(lgd,'location','southeast','fontsize',15);

get(t1,'interpreter');
set(t1,'interpreter','none');
get(t2,'interpreter');
set(t2,'interpreter','none');
% get(t3,'interpreter');
% set(t3,'interpreter','none');

CARDID = param.CARDID;
EIGPATH = param.disperspath;
export_fig(fig_lw,[EIGPATH,CARDID,'.',TYPEID,'.',num2str(N_modes),'_',num2str(Tlim(1)),'_',num2str(Tlim(2)),'s_','lw.pdf'],'-pdf','-painters');
export_fig(fig_gv,[EIGPATH,CARDID,'.',TYPEID,'.',num2str(N_modes),'_',num2str(Tlim(1)),'_',num2str(Tlim(2)),'s_','gv.pdf'],'-pdf','-painters');
export_fig(fig_cv,[EIGPATH,CARDID,'.',TYPEID,'.',num2str(N_modes),'_',num2str(Tlim(1)),'_',num2str(Tlim(2)),'s_','cv.pdf'],'-pdf','-painters');
export_fig(fig_cvgv,[EIGPATH,CARDID,'.',TYPEID,'.',num2str(N_modes),'_',num2str(Tlim(1)),'_',num2str(Tlim(2)),'s_','cvgv.pdf'],'-pdf','-painters');

% print(fig_lw,'-painters','-dpdf','-r400',[EIGPATH,CARDID,'.',TYPEID,'.',num2str(N_modes),'lw.pdf']);
% print(fig_gv,'-painters','-dpdf','-r400',[EIGPATH,CARDID,'.',TYPEID,'.',num2str(N_modes),'gv.pdf']);
% print(fig_cv,'-painters','-dpdf','-r400',[EIGPATH,CARDID,'.',TYPEID,'.',num2str(N_modes),'cv.pdf']);

%% COMPARE FIRST 2 BRANCHES CV and GV
% clr = copper(4);
% fig_gv_v_cv = figure(4); subplot(2,1,2); hold on; box on;
% plot(dat{1}(:,9),dat{1}(:,8),'color',clr(1,:),'linewidth',2);
% plot(dat{2}(:,9),dat{2}(:,8),'color',clr(3,:),'linewidth',2);
% plot(dat{1}(:,10),dat{1}(:,7),'-.','color',clr(1,:),'linewidth',2);
% plot(dat{2}(:,10),dat{2}(:,7),'-.','color',clr(3,:),'linewidth',2);
% legend('Phase T0','Phase T1','Group T0','Group T1','location','southeastoutside');
% axis([5 50 3.5 5.5]);
% xlabel('Period T (sec)','fontsize',12);
% ylabel('Velocity (km/s)','fontsize',12);
% t4 = title(sprintf('%s (maxF = %.1f maxL = %d)',CARDID,maxF,maxL),'fontsize',12);
% get(t4,'interpreter');
% set(t4,'interpreter','none');

%print(fig_gv_v_cv,'-painters','-dpdf','-r400',[EIGPATH,CARDID,'.',TYPEID,'.',num2str(N_modes),'gv_v_cv.pdf']);


