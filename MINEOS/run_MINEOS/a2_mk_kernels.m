%% Driver to Calculate Frechet Kernels in terms of Phase Velocity
% NJA, 2014
% 
% This involves calling fortran programs plot_wk, frechet, frechet_gv,
% frechet_pv
% pylin.patty 2014
%
% JOSH 8/25/2015
%


clear; close all;

parameter_FRECHET;
branch = 0; % Fundamental -> 0

is_deletefrech = 1; % Delete the .frech files to save space?

TYPE = param.TYPE;
CARDID = param.CARDID;

titlename = '';

if ( TYPE == 'T') 
    TYPEID = param.TTYPEID;
elseif ( TYPE == 'S') 
    TYPEID = param.STYPEID;
end

periods = param.periods;

yaxis = [0 350]; %[0 100]; %[0 350];

is_frech_x = 0; % 1 => scale ax; 0 => autoscale
frech_x = [0 2e-8]; %[0 2.0e-7]; %[0 2.0e-8];
%frech_x = [-2.5e-8 2.5e-8];

isfigure = 1; % plot kernels?

%% Set path to executables
setpath_plotwk;

%% Change environment variables to deal with gfortran
setenv('GFORTRAN_STDIN_UNIT', '5') 
setenv('GFORTRAN_STDOUT_UNIT', '6') 
setenv('GFORTRAN_STDERR_UNIT', '0')


%% run plot_wk on the table_hdr file to generate the branch file
write_plotwk(TYPE,CARDID);

com = ['cat run_plotwk.',lower(TYPE),' | plot_wk > plot_wk.LOG'];
[status,log] = system(com);

if status ~= 0     
    error( 'something is wrong at plot_wk')
end


%% run "frechet" to generate the frechet file 
NDISC = 0;
ZDISC = [];

com = ['ls ',param.TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'_1.eig_fix | cat'];
[status eig_fils] = system(com);
if strcmp(eig_fils(end-25:end-1),'No such file or directory')
    disp('Found no *.eig_fix files')
    write_frechet(TYPE,CARDID,NDISC,ZDISC)
else
    disp('Found *.eig_fix files')
    write_frech_chk(NDISC)
end
disp('Be patient! This will take ~25 s');
tic
com = ['cat run_frechet.',lower(TYPE),' | frechet > frechet.LOG'];
[status,log] = system(com);
if status ~= 0     
    error( 'something is wrong at frechet')
end
toc

%% run "frechet_cv" to generate the cv kernels
% Convert frechet to ascii
% Make CV Frechet Kernels
disp('--- Make CV Frechet Kernels ---');


write_frechcv(TYPE,CARDID,branch)

com = ['cat run_frechcv.',lower(TYPE),' | frechet_cv > frechet_cv.LOG'];
[status,log] = system(com);
if status ~= 0     
    error( 'something is wrong at frechet_cv')
end

%% load CARD file (vmod)

CARD = param.CARD;
CARDPATH = param.CARDPATH;
FULLPATH = [CARDPATH,CARD];

fid = fopen(FULLPATH);

%skip 3 line header
for i=1:3
    fgetl(fid);
end

ncard = textscan(fid, '%f%f%f%f%f%f%f%f%f');
fclose(fid);

ncard_temp = ncard;
R = ncard{1};
RHO = ncard{2};
VPV = ncard{3};
VSV = ncard{4};
QKAPPA = ncard{5};
QSHEAR = ncard{6};
VPH = ncard{7};
VSH = ncard{8};
eta = ncard{9};

vs = VSV/1000;
vp = VPV/1000;
r = 6371-R/1000;
%% Convert CV Frechet kernels to ascii with phase-velocity sensitivity
% Will do this for all periods of interest

disp('--- Convert Frechet CV to ascii ---');

    % Program writes run file for draw_frechet_gv, runs it, and reads in
    % sensitivity kernels for all periods of interest
    
if ( TYPE == 'S') 
    FRECH_S = frechcv_asc(TYPE,CARDID,branch);
    if isfigure
        fig1 = figure(62); set(gcf, 'Color', 'w');
        set(gcf,'position',[112   169   830   532]);
        clf
         CC=lines(length(periods));
         CC = flip(brewermap(length(periods),'Spectral'));
         %CC=copper(length(periods));
            ax1 = subplot(1,3,1);
            plot(VSV/1000,r,'linewidth',3,'color',[0 0 0]); hold on;
            plot(VSH/1000,r,'linewidth',3,'color',[1 0 0]);
            set(gca,'Ydir','reverse','linewidth',2,'YMinorTick','on','XMinorTick','on');
            ylim([yaxis]);
            xlim([3 5.2]);
            xlabel('V_{S} (km/s)','fontsize',18);
            ylabel('Depth (km)','fontsize',18);
            legend({'V_{SV}','V_{SH}'},'location','southwest');
            set(gca,'fontsize',16);
%             title('Vs');
            dx = 0.06;
            ax1.Position = [ax1.Position(1:2) ax1.Position(3)+dx ax1.Position(4)];

            for ip = 1:length(periods)
                ax2 = subplot(1,3,2);
                hold on
                dr = gradient(FRECH_S(ip).rad);
                plot(FRECH_S(ip).vsv .* dr,(6371000-FRECH_S(ip).rad)./1000,'-k','linewidth',3,'color',CC(ip,:))
                title('FRECH S - Vsv','fontname','Times New Roman','fontsize',12);
                lgd{ip}=[num2str(periods(ip)),'s'];
                set(gca,'Ydir','reverse','linewidth',2,'YMinorTick','on','XMinorTick','on');
                ylim(yaxis)
                if is_frech_x == 1
                    xlim(frech_x);
                end

                box on;
    %            pause;

            end            
            fig=gcf;
            set(gca,'fontsize',16);
            ax2.Position = [ax1.Position(1)+ax1.Position(3)+0.1 ax2.Position(2) ax2.Position(3)+dx ax2.Position(4)];
            legend(lgd,'position',[ax2.Position(1)+ax2.Position(3)+0.08 0.5 0 0],'box','off');

            
            
            figure(99); clf;
            set(gcf,'position',[119          53        1533         962],'color','w');
           
            subplot(2,3,1); box on; hold on;
            for ip = 1:length(periods)
                plot(FRECH_S(ip).vpv .* dr,(6371000-FRECH_S(ip).rad)./1000,'-k','linewidth',3,'color',CC(ip,:))
            end   
            set(gca,'fontsize',15,'linewidth',1.5,'ydir','reverse');
            xlabel('K_{V_{PV}}')
            ylabel('Depth (km)');
            ylim(yaxis)
            legend(lgd,'location','southeast');
            
            subplot(2,3,2); box on; hold on;
            for ip = 1:length(periods)
                plot(FRECH_S(ip).vph .* dr,(6371000-FRECH_S(ip).rad)./1000,'-k','linewidth',3,'color',CC(ip,:))
            end   
            set(gca,'fontsize',15,'linewidth',1.5,'ydir','reverse');
            xlabel('K_{V_{PH}}')
            ylabel('Depth (km)');
            ylim(yaxis)
            
            subplot(2,3,3); box on; hold on;
            for ip = 1:length(periods)
                plot(FRECH_S(ip).rho .* dr,(6371000-FRECH_S(ip).rad)./1000,'-k','linewidth',3,'color',CC(ip,:))
            end   
            set(gca,'fontsize',15,'linewidth',1.5,'ydir','reverse');
            xlabel('K_{\rho}')
            ylabel('Depth (km)');
            ylim(yaxis)
            
            subplot(2,3,4); box on; hold on;
            for ip = 1:length(periods)
                plot(FRECH_S(ip).vsv .* dr,(6371000-FRECH_S(ip).rad)./1000,'-k','linewidth',3,'color',CC(ip,:))
            end   
            set(gca,'fontsize',15,'linewidth',1.5,'ydir','reverse');
            xlabel('K_{V_{SV}}')
            ylabel('Depth (km)');
            ylim(yaxis)
            
            subplot(2,3,5); box on; hold on;
            for ip = 1:length(periods)
                plot(FRECH_S(ip).vsh .* dr,(6371000-FRECH_S(ip).rad)./1000,'-k','linewidth',3,'color',CC(ip,:))
            end   
            set(gca,'fontsize',15,'linewidth',1.5,'ydir','reverse');
            xlabel('K_{V_{SH}}')
            ylabel('Depth (km)');
            ylim(yaxis)
            
            subplot(2,3,6); box on; hold on;
            for ip = 1:length(periods)
                plot(FRECH_S(ip).eta .* dr,(6371000-FRECH_S(ip).rad)./1000,'-k','linewidth',3,'color',CC(ip,:))
            end   
            set(gca,'fontsize',15,'linewidth',1.5,'ydir','reverse');
            xlabel('K_{\eta}')
            ylabel('Depth (km)');
            ylim(yaxis)

    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ( TYPE == 'T')
    FRECH_T = frechcv_asc(TYPE,CARDID,branch);  

    if isfigure
        fig1 = figure(62); set(gcf, 'Color', 'w'); 
%         set(gcf,'position',[360   165   560   532]);
        set(gcf,'position',[112   169   830   532]);
        clf
%          CC=jet(length(periods));
%          CC=lines(length(periods));
        CC = flip(brewermap(length(periods),'Spectral'));
         %CC=copper(length(periods));
            ax1 = subplot(1,3,1);
            plot(VSV/1000,r,'linewidth',3,'color',[0 0 0]); hold on;
            plot(VSH/1000,r,'linewidth',3,'color',[1 0 0]);
            set(gca,'Ydir','reverse','linewidth',2,'YMinorTick','on','XMinorTick','on');
            ylim([yaxis]);
            xlim([3 5.2]);
            xlabel('V_{S} (km/s)','fontsize',18);
            ylabel('Depth (km)','fontsize',18);
            legend({'V_{SV}','V_{SH}'},'location','southwest');
            set(gca,'fontsize',16);
            dx = 0.06;
            ax1.Position = [ax1.Position(1:2) ax1.Position(3)+dx ax1.Position(4)];
            

            for ip = 1:length(periods)
                ax2 = subplot(1,3,2);
                axis tight;
                hold on
                dr = gradient(FRECH_T(ip).rad);
                plot(FRECH_T(ip).vsh .* dr,(6371000-FRECH_T(ip).rad)./1000,'-k','linewidth',3,'color',CC(ip,:))
%                 title(titlename,'fontsize',18);
                lgd{ip}=[num2str(periods(ip)),'s'];
                set(gca,'YDir','reverse','linewidth',2,'YMinorTick','on','XMinorTick','on')
                ylim(yaxis)
                xlabel('dc/dV_{SH}');
                if is_frech_x == 1
                    xlim(frech_x);
                end

                box on;

            end
            
            fig=gcf;
            set(gca,'fontsize',16);
            ax2.Position = [ax1.Position(1)+ax1.Position(3)+0.1 ax2.Position(2) ax2.Position(3)+dx ax2.Position(4)];
            legend(lgd,'position',[ax2.Position(1)+ax2.Position(3)+0.08 0.5 0 0],'box','off');

            
            figure(99); clf;
            set(gcf,'position',[119          53        1533         962/2],'color','w');
            
            subplot(1,3,1); box on; hold on;
            for ip = 1:length(periods)
                plot(FRECH_T(ip).vsv .* dr,(6371000-FRECH_T(ip).rad)./1000,'-k','linewidth',3,'color',CC(ip,:))
            end   
            set(gca,'fontsize',15,'linewidth',1.5,'ydir','reverse');
            xlabel('K_{V_{SV}}')
            ylabel('Depth (km)');
            ylim(yaxis)
            
            subplot(1,3,2); box on; hold on;
            for ip = 1:length(periods)
                plot(FRECH_T(ip).vsh .* dr,(6371000-FRECH_T(ip).rad)./1000,'-k','linewidth',3,'color',CC(ip,:))
            end   
            set(gca,'fontsize',15,'linewidth',1.5,'ydir','reverse');
            xlabel('K_{V_{SH}}')
            ylabel('Depth (km)');
            ylim(yaxis)
            
            subplot(1,3,3); box on; hold on;
            for ip = 1:length(periods)
                plot(FRECH_T(ip).rho .* dr,(6371000-FRECH_T(ip).rad)./1000,'-k','linewidth',3,'color',CC(ip,:))
            end   
            set(gca,'fontsize',15,'linewidth',1.5,'ydir','reverse');
            xlabel('K_{\rho}')
            ylabel('Depth (km)');
            ylim(yaxis)
            
    end
end

FRECHETPATH = param.frechetpath;
delete(['run_plotwk.',lower(TYPE)],['run_frechcv.',lower(TYPE)],['run_frechet.',lower(TYPE)],['run_frechcv_asc.',lower(TYPE)]);
save2pdf([FRECHETPATH,'CARD_Vs_kernels_',lower(TYPE),'_',CARDID,'_b',num2str(branch),'.',num2str(N_modes),'_',num2str(periods(1)),'_',num2str(periods(end)),'s.pdf'],fig1,1000)

%    savefile = [CARD,'_fcv.mat'];
%    save(savefile,'FRECH_T','FRECH_S');
% Change the environment variables back to the way they were
setenv('GFORTRAN_STDIN_UNIT', '-1') 
setenv('GFORTRAN_STDOUT_UNIT', '-1') 
setenv('GFORTRAN_STDERR_UNIT', '-1')

delete('*.LOG');
if is_deletefrech
    delete([param.frechetpath,'*.fcv.*'],[param.frechetpath,'*.frech'])
end