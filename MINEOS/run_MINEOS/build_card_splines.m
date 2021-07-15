%%% 7/6/16 -- JBR
%
% Program to build card splines
% Takes 6 knots in the upper 300 km of starting model and varies the top 5 knots from 0% to 3%
% anisotropy, keeping the deepest one fixed. Therefore there are 2^5 = 32
% total models generated.
%
% COLUMNS: 1   2  3   4     5     6     7   8   9
%          R,RHO,VPV,VSV,QKAPPA,QSHEAR,VPH,VSH,ETA
%


clear all; close all;
parameter_FRECHET

CARDID = param.CARDID;
frechetpath = param.frechetpath;
figname = [CARDID,'_3perc'];

YLIMS = [0 400];

%knot = [760   750   740   730   720   710   703];
knot = [760   750   740   730   720    703];
aniso_strength = 1.03; % 3 percent
eta_aniso = 0.95;

FIGPATH = [param.CARDPATH,'/SPLINECARDS/'];
system(['mkdir ',FIGPATH]);
FIGPATH = [FIGPATH,CARDID,'/'];
system(['mkdir ',FIGPATH]);


%%% MAKE DIRECTORY FOR SPLINE CARDS %%%
SPLINEPATH = [param.CARDPATH,'/',param.CARDID,'/'];
system(['mkdir ',SPLINEPATH]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_ncard = [param.CARDPATH,param.CARD];
h_nqmod = [param.CARDPATH,param.CARDID,'.qmod'];

fidr = fopen(h_ncard,'r');
hdr_l1 = fgetl(fidr);
hdr_l2 = fgetl(fidr);
hdr_l3 = fgetl(fidr);
ncard = textscan(fidr, '%f%f%f%f%f%f%f%f%f');
fclose(fidr);

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

% CALCULATE ANISOTROPY AT EACH KNOT
VSH_aniso = VSH(knot)*aniso_strength;
VPH_aniso = VPH(knot)*aniso_strength;
eta_aniso = ones(size(VSH))*eta_aniso;


% LOOP THROUGH THE MODEL SPACE
% fig1 = figure(1); clf; hold on; 
% set(gcf, 'Color', 'w');
% set(gcf,'position',[360   165   560   532]);
% fig2 = figure(2); clf; hold on; 
% set(gcf, 'Color', 'w');
% set(gcf,'position',[360   165   560   532]);

nknots = length(knot);
nloops = 2^(nknots-1); % keep deepest knot fixed
clr = jet(nloops);

iicard = 0;
%test = [];
for iknot = 1:nknots
    if iknot == 1 % set first model to starting model
        I_knot_aniso = logical(zeros(1,nknots-1));
        I_aniso = I_knot_aniso;
        VSH_new = VSH(knot);
        VSH_new(I_aniso) = VSH_aniso(I_aniso);
        VPH_new = VPH(knot);
        VPH_new(I_aniso) = VPH_aniso(I_aniso);
        eta_new = eta(knot);
        eta_new(I_aniso) = eta_aniso(I_aniso);
        % CALCULATE SPLINES
        y = R(knot);
        x_VSH = VSH_new;
        x_VPH = VPH_new;
        x_eta = eta_new;
        yy = R(knot(end):knot(1));
        VSH_spline = spline(y,x_VSH,yy);
        VSH_pchip = pchip(y,x_VSH,yy); % Piecewise Cubic Hermite Interpolating Polynomial
        VPH_pchip = pchip(y,x_VPH,yy);
        eta_pchip = pchip(y,x_eta,yy);
        iicard = iicard + 1;
        
%         % PLOT STARTING MODEL
%         figure(1);
%         subplot(1,2,1); hold on;
%         plot(VSH/1000,6371-R/1000,'k','linewidth',2)
%         plot(VSH(knot)/1000,6371-R(knot)/1000,'ok','linewidth',2)
%         
%         % PLOT SPLINE
%         plot(VSH_pchip/1000,6371-yy/1000,'-','color',clr(iicard,:));
%         
%         % PLOT SINGLE MODEL
%         figure(2); clf;
%         subplot(1,2,1); hold on; box on;
%         plot(VSH/1000,6371-R/1000,'-k','linewidth',2)
%         plot(VSH(knot)/1000,6371-R(knot)/1000,'ok','linewidth',2)
%         
%         plot(VSH_pchip/1000,6371-yy/1000,'-r','linewidth',2);
%         set(gca,'Ydir','reverse');
%         ylim(YLIMS);
%         xlim([4 5]);
%         %ylim([4 7]);
%         %xlim([0 4])
%         xlabel('V_{SH} (km/s)','fontsize',18);
%         ylabel('Depth (km)','fontsize',18);
%         set(gca,'fontsize',18);
%         
%         outfil_name = [FIGPATH,figname,'_',num2str(iicard),'.pdf'];
%         save2pdf([outfil_name],fig2,1000);

        VSH_out = VSH;
        VSH_out(knot(end):knot(1)) = VSH_pchip;
        VPH_out = VPH;
        VPH_out(knot(end):knot(1)) = VPH_pchip;
        eta_out = eta;
        eta_out(knot(end):knot(1)) = eta_pchip;
        
        ocard{1} = ncard{1}; % R
        ocard{2} = ncard{2}; % RHO
        ocard{3} = VPV; % VPV
        ocard{4} = VSV; % VSV
        ocard{5} = ncard{5}; % QKAPPA
        ocard{6} = ncard{6}; % QSHEAR 
        ocard{7} = VPH_out; % VPH 
        ocard{8} = VSH_out; % VSH
        ocard{9} = eta_out; % ETA

        % specify number of rows in model
        numrow = size(ocard{1},1);
        hdr_l3 = [' ',num2str(numrow),' ', hdr_l3(6:end)];
        
        OCARD = [param.CARDID,'_',num2str(iicard)];
        h_ocard = [SPLINEPATH,OCARD,'.card'];
        fidw = fopen(h_ocard,'w');
        fprintf(fidw,'%s\n',hdr_l1);
        fprintf(fidw,'%s\n',hdr_l2);
        fprintf(fidw,'%s\n',hdr_l3);
        for i = 1:numrow
            fprintf(fidw,'%7.0f.%9.2f%9.2f%9.2f%9.1f%9.1f%9.2f%9.2f%9.5f\n',ocard{1}(i),ocard{2}(i),ocard{3}(i),ocard{4}(i),ocard{5}(i),ocard{6}(i),ocard{7}(i),ocard{8}(i),ocard{9}(i));
        end
        fclose(fidw);

        h_oqmod = [SPLINEPATH,OCARD,'.qmod'];
        com = ['cp ',h_nqmod,' ',h_oqmod];
        system(com);

    else
        I_knot_aniso(iknot-1) = true;
        I_mat_aniso = unique(perms(I_knot_aniso),'rows');
        nperms = size(I_mat_aniso,1);
        for iperm = 1:nperms
            I_aniso = logical(I_mat_aniso(iperm,:));
            %test = [test; I_aniso];
            VSH_new = VSH(knot);
            VSH_new(I_aniso) = VSH_aniso(I_aniso);
            VPH_new = VPH(knot);
            VPH_new(I_aniso) = VPH_aniso(I_aniso);
            eta_new = eta(knot);
            eta_new(I_aniso) = eta_aniso(I_aniso);
            % CALCULATE SPLINES
            y = R(knot);
            x_VSH = VSH_new;
            x_VPH = VPH_new;
            x_eta = eta_new;
            yy = R(knot(end):knot(1));
            VSH_spline = spline(y,x_VSH,yy);
            VSH_pchip = pchip(y,x_VSH,yy); % Piecewise Cubic Hermite Interpolating Polynomial
            VPH_pchip = pchip(y,x_VPH,yy);
            eta_pchip = pchip(y,x_eta,yy);
            iicard = iicard + 1;
            
%             % PLOT SPLINE
%             figure(1);
%             subplot(1,2,1);
%             plot(VSH_pchip/1000,6371-yy/1000,'-','color',clr(iicard,:));
%             
%             %pause;
%             
%             % PLOT SINGLE MODEL
%             figure(2); clf;
%             subplot(1,2,1); hold on; box on;
%             plot(VSH/1000,6371-R/1000,'-k','linewidth',2)
%             plot(VSH(knot)/1000,6371-R(knot)/1000,'ok','linewidth',2)
% 
%             plot(VSH_pchip/1000,6371-yy/1000,'-r','linewidth',2);
%             set(gca,'Ydir','reverse');
%             ylim(YLIMS);
%             xlim([4 5]);
%             %ylim([4 7]);
%             %xlim([0 4])
%             xlabel('V_{SH} (km/s)','fontsize',18);
%             ylabel('Depth (km)','fontsize',18);
%             set(gca,'fontsize',18);
% 
%             outfil_name = [FIGPATH,figname,'_',num2str(iicard),'.pdf'];
%             save2pdf([outfil_name],fig2,1000);

            VSH_out = VSH;
            VSH_out(knot(end):knot(1)) = VSH_pchip;
            VPH_out = VPH;
            VPH_out(knot(end):knot(1)) = VPH_pchip;
            eta_out = eta;
            eta_out(knot(end):knot(1)) = eta_pchip;

            ocard{1} = ncard{1}; % R
            ocard{2} = ncard{2}; % RHO
            ocard{3} = VPV; % VPV
            ocard{4} = VSV; % VSV
            ocard{5} = ncard{5}; % QKAPPA
            ocard{6} = ncard{6}; % QSHEAR 
            ocard{7} = VPH_out; % VPH 
            ocard{8} = VSH_out; % VSH
            ocard{9} = eta_out; % ETA

            % specify number of rows in model
            numrow = size(ocard{1},1);
            hdr_l3 = [' ',num2str(numrow),' ', hdr_l3(6:end)];
            
            OCARD = [param.CARDID,'_',num2str(iicard)];
            h_ocard = [SPLINEPATH,OCARD,'.card'];
            fidw = fopen(h_ocard,'w');
            fprintf(fidw,'%s\n',hdr_l1);
            fprintf(fidw,'%s\n',hdr_l2);
            fprintf(fidw,'%s\n',hdr_l3);
            for i = 1:numrow
                fprintf(fidw,'%7.0f.%9.2f%9.2f%9.2f%9.1f%9.1f%9.2f%9.2f%9.5f\n',ocard{1}(i),ocard{2}(i),ocard{3}(i),ocard{4}(i),ocard{5}(i),ocard{6}(i),ocard{7}(i),ocard{8}(i),ocard{9}(i));
            end
            fclose(fidw);

            h_oqmod = [SPLINEPATH,OCARD,'.qmod'];
            com = ['cp ',h_nqmod,' ',h_oqmod];
            system(com);
        end
    end
end
% figure(1);
% subplot(1,2,1); box on;
% set(gca,'Ydir','reverse');
% ylim(YLIMS);
% xlim([4 5]);
% %ylim([4 7]);
% %xlim([0 4])
% xlabel('V_{SH} (km/s)','fontsize',18);
% ylabel('Depth (km)','fontsize',18);
% set(gca,'fontsize',18);
% 
% outfil_name = [FIGPATH,figname,'_allcards.pdf'];
% save2pdf([outfil_name],fig1,1000);
