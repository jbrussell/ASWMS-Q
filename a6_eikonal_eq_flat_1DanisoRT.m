% Read in the eventcs structures and apply eikonal tomography on each event.
% Written by Ge Jin, jinwar@gmail.com
% 2013.1.16
%
% JBR 3/19/19 : add first derivative smoothing 
% JBR 11/19 : Invert all events directly for isotropic velocity and
% 1-D azimuthal anisotropy.
%
clear
%plot native

% setup parameters
setup_parameters
setup_ErrorCode

% JBR
cohere_tol = parameters.cohere_tol;

% % Smoothing parameters
% flweight_array = 100*ones(length(parameters.periods)); %parameters.flweight_array
% dterrtol = 2;    % largest variance of the inversion error allowed
% inverse_err_tol = 2; %2  % count be number of standard devition
% azi_bin_deg = 20; % [deg] size of azimuthal bins
% min_nbin = 10; % minimum number of measurements in order to include bin
% % Main QC parameters
% fiterr_tol = 1e-2; % wavelet fitting error, throw out measurements greater than this
% maxstadist = 600;
% minstadist = 200;
% cohere_tol = 0.80; % 0.65
% min_stadist_wavelength = 0.33; %0.5; % minimum station separation in wavelengths
% max_stadist_wavelength = 999;
% ref_phv = [3.9936 4.0041 4.0005 3.9999 3.9929 3.9832 3.9813 3.9841 3.9874 3.9996 4.0138 4.0519 4.0930 4.1677 4.2520]; % for calculating wavelength
% APM = 114; % absolute plate motion (GSRM 2.1; NNR) https://www.unavco.org/software/geodetic-utilities/plate-motion-calculator/plate-motion-calculator.html
% FSD = 75; 

% Smoothing parameters
flweight_array = 100*ones(length(parameters.periods)); %parameters.flweight_array
dterrtol = 2;    % largest variance of the inversion error allowed
inverse_err_tol = 2; %2  % count be number of standard devition
azi_bin_deg = 20; % [deg] size of azimuthal bins
min_nbin = 10; % minimum number of measurements in order to include bin
% Main QC parameters
fiterr_tol = 1e-2; % wavelet fitting error, throw out measurements greater than this
maxstadist = 600;
minstadist = 200;
cohere_tol = 0.80; % 0.65
min_stadist_wavelength = 0.33; %0.5; % minimum station separation in wavelengths
max_stadist_wavelength = 999;
ref_phv = [3.9936 4.0041 4.0005 3.9999 3.9929 3.9832 3.9813 3.9841 3.9874 3.9996 4.0138 4.0519 4.0930 4.1677 4.2520]; % for calculating wavelength
APM = 114; % absolute plate motion (GSRM 2.1; NNR) https://www.unavco.org/software/geodetic-utilities/plate-motion-calculator/plate-motion-calculator.html
FSD = 75; 

% Norm damping for azimuthal anisotropy
damp_azi = [1 1 1e10 1e10]; % [2c 2s 4c 4s] % Damping individual parameters
aziweight = 1; % global weight

% debug setting
isfigure = 1;
isdisp = 1;
is_overwrite = 1;

% % input path
% eventcs_path = './CSmeasure/';
% % output path
% eikonl_ani_output_path = './eikonal/';

workingdir = parameters.workingdir;
% input path
eventcs_path = [workingdir,'CSmeasure/'];
% output path
eikonl_output_path = [workingdir,'eikonal/'];
eikonl_ani_output_path = [workingdir];


comp = parameters.component;
lalim=parameters.lalim;
lolim=parameters.lolim;
gridsize=parameters.gridsize;
periods = parameters.periods;
raydensetol=parameters.raydensetol;
smweight_array = parameters.smweight_array;
% flweight_array = parameters.flweight_array; % JBR
Tdumpweight0 = parameters.Tdumpweight;
Rdumpweight0 = parameters.Rdumpweight;
fiterrtol = parameters.fiterrtol;
% dterrtol = parameters.dterrtol;
isRsmooth = parameters.isRsmooth;
% inverse_err_tol = parameters.inverse_err_tol;
min_amp_tol  = parameters.min_amp_tol;

% setup useful variables
xnode=lalim(1):gridsize:lalim(2);
ynode=lolim(1):gridsize:lolim(2);
Nx=length(xnode);
Ny=length(ynode);
[xi yi]=ndgrid(xnode,ynode);

% Setup universal smoothing kernel
disp('initial the smoothing kernel')
tic
	% longtitude smoothing
    [i,j] = ndgrid(1:Nx,2:(Ny-1));
    ind = j(:) + Ny*(i(:)-1);
    dy = diff(ynode)*cosd(mean(xnode));  % correct smoothing for latitude
    dy1 = dy(j(:)-1);
    dy2 = dy(j(:));

    Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
                    [-2./(dy1.*(dy1+dy2)), 2./(dy1.*dy2), -2./(dy2.*(dy1+dy2))],Nx*Ny+4,Nx*Ny+4);

	% latitude smoothing
    [i,j] = ndgrid(2:(Nx-1),1:Ny);
    ind = j(:) + Ny*(i(:)-1);
    dx = diff(xnode);
    dx1 = dx(i(:)-1);
    dx2 = dx(i(:));

    Areg = [Areg;sparse(repmat(ind,1,3),[ind-Ny,ind,ind+Ny], ...
            [-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), -2./(dx2.*(dx1+dx2))],Nx*Ny+4,Nx*Ny+4)];

    F=sparse(Nx*Ny*2*2,Nx*Ny*2+4);
    for n=1:size(Areg,1)
        ind=find(Areg(n,:)~=0);
        F(2*n-1,2*ind-1)=Areg(n,ind);
        F(2*n,2*ind)=Areg(n,ind);
    end
toc

% JBR - define first derivative "flatness" kernel
F2 = flat_kernel_build(xnode, ynode, Nx*Ny);
Ftemp = sparse(size(F2,1),size(F2,2)+4);
Ftemp(1:Nx*Ny*2*2,1:Nx*Ny*2) = F2;
F2 = Ftemp;

%JBR - add anisotropic terms to smoothing matrix (add 4 columns and 4 rows)
% damping matrix
Areg_azi = zeros(4,Nx*Ny*2+4);
azipart = eye(4,4) .* diag(damp_azi);
Areg_azi(1:4,Nx*Ny*2+1:Nx*Ny*2+4) = azipart;
F_azi_damp = Areg_azi;

% read in bad station list, if existed
if exist('badsta.lst')
	badstnms = textread('badsta.lst','%s');
	disp('Found Bad stations:')
	disp(badstnms)
end

csmatfiles = dir([eventcs_path,'/*cs_',comp,'.mat']);
% Loop through the periods
clear eventphv_ani;
for ip = 1:length(periods)
	clear mat_azi azi rays ddist dt w 
    disp(periods(ip));
	smweight0 = smweight_array(ip);
	flweight0 = flweight_array(ip); % JBR
 	raynum = 0;
	for ie = 1:length(csmatfiles)
	%for ie = 30
		% read in data and set up useful variables
		temp = load([eventcs_path,csmatfiles(ie).name]);
		eventcs =  temp.eventcs;
% 		disp(eventcs.id)
		evla = eventcs.evla;
		evlo = eventcs.evlo;

		matfilename = [eikonl_ani_output_path,'/eikonal_ani_',comp,'.mat'];
		if exist(matfilename,'file') && ~is_overwrite
			disp(['Exist ',matfilename,', skip']);
			continue;
        end
        
%        % Load wavefield azimuth from previous eikonal solution
%         matfilename_in = [eikonl_output_path,'/',eventcs.id,'_eikonal_',comp,'.mat'];
%         if exist(matfilename_in,'file')==2
%             temp = load(matfilename_in);
%             azi_ev = angle( nanmedian(temp.eventphv(ip).GVx(:)) + nanmedian(temp.eventphv(ip).GVy(:)).*sqrt(-1));
%             azi_ev = rad2deg(azi_ev);
% %             azi_ev = nanmedian(azi_ev(:));
%             azi_flag = 1;
%             if isnan(azi_ev)
%                 azi_flag = 0;
%             end
%         else
%             azi_flag = 0;
%         end


		if exist('badstnms','var')
			badstaids = find(ismember(eventcs.stnms,badstnms));
		else
			badstaids = [];
		end

		% Build the rotation matrix
		razi = azimuth(xi+gridsize/2,yi+gridsize/2,evla,evlo)+180;
		R = sparse(2*Nx*Ny+4,2*Nx*Ny+4);
		for i=1:Nx
			for j=1:Ny
				n=Ny*(i-1)+j;
				theta = razi(i,j);
				R(2*n-1,2*n-1) = cosd(theta);
				R(2*n-1,2*n) = sind(theta);
				R(2*n,2*n-1) = -sind(theta);
				R(2*n,2*n) = cosd(theta);
			end
		end

		% Calculate the relative travel time compare to one reference station
% 		travel_time = Cal_Relative_dtp(eventcs);

		% Build the ray locations
		for ics = 1:length(eventcs.CS)
			raynum = raynum+1;
			
			if eventcs.CS(ics).cohere(ip)<cohere_tol && eventcs.CS(ics).isgood(ip)>0
                eventcs.CS(ics).isgood(ip) = ErrorCode.low_cohere;
            end
			if (eventcs.CS(ics).ddist < ref_phv(ip)*periods(ip)*min_stadist_wavelength || ...
                    eventcs.CS(ics).ddist > ref_phv(ip)*periods(ip)*max_stadist_wavelength) && eventcs.CS(ics).isgood(ip)>0
				eventcs.CS(ics).isgood(ip) = ErrorCode.min_stadist_wavelength;
            end
            if (eventcs.CS(ics).ddist > maxstadist || ...
                    eventcs.CS(ics).ddist < minstadist) && eventcs.CS(ics).isgood(ip)>0
                eventcs.CS(ics).isgood(ip) = -13;
            end
            if (eventcs.CS(ics).fiterr(ip) > fiterr_tol)
                eventcs.CS(ics).isgood(ip) = -14;
            end
            if eventcs.CS(ics).cohere(ip)>=cohere_tol && eventcs.CS(ics).isgood(ip)==ErrorCode.low_cohere
                eventcs.CS(ics).isgood(ip) = 1;
            end
			if eventcs.CS(ics).isgood(ip) > 0 
				dt(raynum,:) = eventcs.CS(ics).dtp(ip);
				w(raynum,:) = 1;
%                 w(raynum,:) = eventcs.CS(ics).fiterr(ip)^(-1/2);
%                 w(raynum,:) = eventcs.CS(ics).fiterr(ip)^(-1/2) * eventcs.CS(ics).ddist^(1/2);
			else
				dt(raynum,:) = eventcs.CS(ics).dtp(ip);
				w(raynum,:) = 0;
			end
			if sum(ismember([eventcs.CS(ics).sta1 eventcs.CS(ics).sta2],badstaids)) > 0
				w(raynum,:) = 0;
            end
            
            evinfo.evid{raynum,:} = eventcs.id;
            evinfo.evla(raynum,:) = eventcs.evla;
            evinfo.evlo(raynum,:) = eventcs.evlo;
            evinfo.evdp(raynum,:) = eventcs.evdp;
            evinfo.Mw(raynum,:) = eventcs.Mw;
			
			rays(raynum,1) = eventcs.stlas(eventcs.CS(ics).sta1);
			rays(raynum,2) = eventcs.stlos(eventcs.CS(ics).sta1);
			rays(raynum,3) = eventcs.stlas(eventcs.CS(ics).sta2);
			rays(raynum,4) = eventcs.stlos(eventcs.CS(ics).sta2);
		
			% JBR - Build azimuthal part of data kernel
			ddist(raynum,:) = eventcs.CS(ics).ddist;
	        % dt(ics) = eventcs.dtp(ics);
	        % phv(ics) = ddist(ics)./dt(ics);
            mean_stala = mean([eventcs.stlas(eventcs.CS(ics).sta1), eventcs.stlas(eventcs.CS(ics).sta2)]);
            mean_stalo = mean([eventcs.stlos(eventcs.CS(ics).sta1), eventcs.stlos(eventcs.CS(ics).sta2)]);
            
%           % Build azimuth vector. Use wavefield if available, if not use back azimuth from station-pair midpoint
% 	        if azi_flag==1
%                 azi(raynum,:) = azi_ev;
%             else
%                 [~,azi(raynum,:)]=distance(mean_stala,mean_stalo,...
%                                             evla,evlo);
% %                 w(raynum,:) = 0;
%             end
            azi(raynum,:)=azimuth(mean_stala,mean_stalo,...
                                       evla,evlo);
            
% 	        if azi(raynum) > 180
% 	            azi(raynum,:) = azi(raynum) - 360;
% 	        end
	    	mat_azi(raynum,:) = ddist(raynum) * [cosd(2*azi(raynum)), sind(2*azi(raynum)), cosd(4*azi(raynum)), sind(4*azi(raynum)) ];
        end
    end
		W = sparse(diag(w));
		
		% Build the kernel
		disp('Building up ray path kernel')
		tic
			% mat_iso=kernel_build(rays,xnode,ynode);
            mat_iso=kernel_build_RT(rays,xnode,ynode,azi);
		toc
		mat = [mat_iso, mat_azi];

		% build dumping matrix for ST
		dumpmatT = R(2:2:2*Nx*Ny,:);
		
		% build dumping matrix for SR
		dumpmatR = R(1:2:2*Nx*Ny-1,:);
	

		% Normalize smoothing kernel
        NR=norm(F,1);
        NA=norm(W*mat,1);
        smweight = smweight0*NA/NR;
        
        % JBR - Normalize flatness kernel
        NR=norm(F2,1);
        NA=norm(W*mat,1);
        flweight = flweight0*NA/NR;

		% Normalize dumping matrix for ST
		NR=norm(dumpmatT,1);
		NA=norm(W*mat,1);
		dumpweightT = Tdumpweight0*NA/NR;
		
		% Normalize dumping matrix for SR
		NR=norm(dumpmatR,1);
		NA=norm(W*mat,1);
		dumpweightR = Rdumpweight0*NA/NR;
		
		% Rescale azimuthal anisotropy damping
		NR=norm(F_azi_damp,1);
		NA=norm(W*mat,1);
		aziweight0 = aziweight*NA/NR;

		% Set up matrix on both side
		if isRsmooth
			% A=[W*mat;smweight*F*R;flweight*F2*R;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
			A=[W*mat;smweight*F*R;flweight*F2*R;dumpweightT*dumpmatT;dumpweightR*dumpmatR;aziweight0*F_azi_damp];
			% A=[W*mat; smweight*F; aziweight*F_azi_damp];
		else
			% A=[W*mat;smweight*F;flweight*F2;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
			A=[W*mat;smweight*F;flweight*F2;dumpweightT*dumpmatT;dumpweightR*dumpmatR;aziweight0*F_azi_damp];
		end

		avgv = eventcs.avgphv(ip);
        % rhs=[W*dt;zeros(size(F,1),1);zeros(size(F2,1),1);zeros(size(dumpmatT,1),1);dumpweightR*ones(size(dumpmatR,1),1)./avgv];
		rhs=[W*dt;zeros(size(F,1),1);zeros(size(F2,1),1);zeros(size(dumpmatT,1),1);dumpweightR*ones(size(dumpmatR,1),1)./avgv;zeros(size(F_azi_damp,1),1)];

		
		% Least square inversion
        phaseg=(A'*A)\(A'*rhs);
	        
        % Iteratively down weight the measurement with high error
		niter=0;
		ind = find(diag(W)==0);
		if isdisp
			disp(['Before iteration'])
			disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
			disp(['Bad Measurement Number: ', num2str(length(ind))]);
		end
        niter=1;
        while niter < 2
            niter=niter+1;
            isgood_mat = diag(diag(W)>0);
            err = mat*phaseg - dt;
% 			err = W*err;
            err = isgood_mat*err;
%             err = W*err;
            stderr=std(err);
            if stderr > dterrtol
                stderr = dterrtol;
            end
            for i=1:length(err)
                if abs(err(i)) > inverse_err_tol*stderr
                    W(i,i)=0;
                end
            end
            ind = find(diag(W)==0);
			if isdisp
				disp('After:')
				disp(['Good Measurement Number: ', num2str(length(diag(W))-length(ind))]);
				disp(['Bad Measurement Number: ', num2str(length(ind))]);
			end
            
            % Rescale the smooth kernel
            NR=norm(F,1);
            NA=norm(W*mat,1);
            smweight = smweight0*NA/NR;
            
            % JBR - Normalize flatness kernel
            NR=norm(F2,1);
            NA=norm(W*mat,1);
            flweight = flweight0*NA/NR;
            
            % rescale dumping matrix for St
            NR=norm(dumpmatT,1);
            NA=norm(W*mat,1);
            dumpweightT = Tdumpweight0*NA/NR;
            
            % rescale dumping matrix for SR
            NR=norm(dumpmatR,1);
            NA=norm(W*mat,1);
            dumpweightR = Rdumpweight0*NA/NR;
			
			% Rescale azimuthal anisotropy damping
            NR=norm(F_azi_damp,1);
            NA=norm(W*mat,1);
            aziweight0 = aziweight*NA/NR;
            
            if isRsmooth
                % A=[W*mat;smweight*F*R;flweight*F2*R;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
				A=[W*mat;smweight*F*R;flweight*F2*R;dumpweightT*dumpmatT;dumpweightR*dumpmatR;aziweight0*F_azi_damp];
				% A=[W*mat; smweight*F; aziweight*F_azi_damp];
            else
                % A=[W*mat;smweight*F;flweight*F2;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
				A=[W*mat;smweight*F;flweight*F2;dumpweightT*dumpmatT;dumpweightR*dumpmatR;aziweight0*F_azi_damp];
            end
            % rhs=[W*dt;zeros(size(F,1),1);zeros(size(F2,1),1);zeros(size(dumpmatT,1),1);dumpweightR*ones(size(dumpmatR,1),1)./avgv];
			rhs=[W*dt;zeros(size(F,1),1);zeros(size(F2,1),1);zeros(size(dumpmatT,1),1);dumpweightR*ones(size(dumpmatR,1),1)./avgv;zeros(size(F_azi_damp,1),1)];
            phaseg=(A'*A)\(A'*rhs);
        end	

        % Calculate the kernel density
        %sumG=sum(abs(mat),1);
        ind=1:Nx*Ny;
        rayW = W;
        rayW(find(rayW>1))=1;
        raymat = rayW*mat;
        sumG(ind)=sum((raymat(:,2*ind).^2+raymat(:,2*ind-1).^2).^.5,1);
        clear raydense
        for i=1:Nx
            for j=1:Ny
                n=Ny*(i-1)+j;
                raydense(i,j)=sumG(n);
            end
        end
        
        %        disp(' Get rid of uncertainty area');
        fullphaseg = phaseg;
        for i=1:Nx
            for j=1:Ny
                n=Ny*(i-1)+j;
                if raydense(i,j) < raydensetol %&& ~issyntest
                    phaseg(2*n-1)=NaN;
                    phaseg(2*n)=NaN;
                end
            end
        end
        
        % Isotropic phase velocity
        phv_iso = ddist./(mat_iso*phaseg(1:Nx*Ny*2));
        
		% Change phaseg into phase velocity
		for i=1:Nx
			for j=1:Ny
				n=Ny*(i-1)+j;
				GVx(i,j)= phaseg(2*n-1);
				GVy(i,j)= phaseg(2*n);
			end
		end
		GV=(GVx.^2+GVy.^2).^-.5;
        GV_azi = angle(GVx + GVy.*sqrt(-1));
		GV_azi = rad2deg(GV_azi);
		GV_av = 1./nanmean(1./GV(:));
		slow_av = 1./GV_av;
		% get azimuthal anisotropy
		Ac2 = -phaseg(Nx*Ny*2+1)./slow_av;
	    As2 = -phaseg(Nx*Ny*2+2)./slow_av;
	    Ac4 = -phaseg(Nx*Ny*2+3)./slow_av;
	    As4 = -phaseg(Nx*Ny*2+4)./slow_av;
	    A2 = sqrt(Ac2^2+As2^2);
	    A4 = sqrt(Ac4^2+As4^2);
	    phi2 = 1/2*atan2d(As2,Ac2);
	    phi4 = 1/4*atan2d(As4,Ac4);

		% save the result in the structure
		eventphv_ani(ip).rays = rays;
		eventphv_ani(ip).w = diag(W);
		eventphv_ani(ip).goodnum = length(find(eventphv_ani(ip).w>0));
		eventphv_ani(ip).badnum = length(find(eventphv_ani(ip).w==0));
		eventphv_ani(ip).dt = dt;
		eventphv_ani(ip).GV = GV;
        eventphv_ani(ip).GV_azi = GV_azi;
        eventphv_ani(ip).GV_av = GV_av;
		eventphv_ani(ip).GVx = GVx;
		eventphv_ani(ip).GVy = GVy;
        eventphv_ani(ip).phv_iso = phv_iso;
		eventphv_ani(ip).raydense = raydense;
		eventphv_ani(ip).lalim = lalim;
		eventphv_ani(ip).lolim = lolim;
		eventphv_ani(ip).gridsize = gridsize;
		eventphv_ani(ip).azi = azi;
		eventphv_ani(ip).rays = rays; 
		eventphv_ani(ip).ddist = ddist;
		eventphv_ani(ip).id = evinfo.evid;
		eventphv_ani(ip).evla = evinfo.evla;
		eventphv_ani(ip).evlo = evinfo.evlo;
		eventphv_ani(ip).evdp = evinfo.evdp;
        eventphv_ani(ip).Mw = evinfo.Mw;
		eventphv_ani(ip).period = periods(ip);
		% eventphv_ani(ip).traveltime = travel_time(ip).tp;
		% eventphv_ani(ip).stlas = eventcs.stlas;
		% eventphv_ani(ip).stlos = eventcs.stlos;
		% eventphv_ani(ip).stnms = eventcs.stnms;
        eventphv_ani(ip).isgood = eventphv_ani(ip).w>0;
        eventphv_ani(ip).phv = ddist./dt;
		eventphv_ani(ip).A2 = A2;
		eventphv_ani(ip).A4 = A4;
		eventphv_ani(ip).phi2 = phi2;
		eventphv_ani(ip).phi4 = phi4;
        
		disp(['Period:',num2str(periods(ip)),', Goodnum:',num2str(eventphv_ani(ip).goodnum),...
				'Badnum:',num2str(eventphv_ani(ip).badnum)]);
	end % end of periods loop
    
%% Make event list
fid = fopen([workingdir,'azi_evlist.txt'],'w');

for ip = 7 % 50 s
    isgood = eventphv_ani(ip).isgood;
    evids_all = eventphv_ani(ip).id(isgood);
    [evids,I] = unique(evids_all);
    evlas = eventphv_ani(ip).evla(isgood); evlas = evlas(I);
    evlos = eventphv_ani(ip).evlo(isgood); evlos = evlos(I);
    evdps = eventphv_ani(ip).evdp(isgood); evdps = evdps(I);
    evMws = eventphv_ani(ip).Mw(isgood);   evMws = evMws(I);
    for iev = 1:length(evids)
        evid = evids{iev};
        evla = evlas(iev);
        evlo = evlos(iev);
        evdp = evdps(iev);
        evMw = evMws(iev);
        fprintf(fid,'%13s %12f %12f %5.1f %5.2f\n',evid,evla,evlo,evdp,evMw);
    end
%     display(length(evids))
end
fclose(fid);
    
 %% Fit anisotropy
fit_azi = [];
for ip = 1:length(periods)
    isgood = eventphv_ani(ip).isgood;
    phv = eventphv_ani(ip).phv(isgood);
    azi = eventphv_ani(ip).azi(isgood);
    w = eventphv_ani(ip).w(isgood);
    w(w > (median(w)+rms(w)) ) = median(w)+rms(w);
    weight = w.*(-1/2);
    [para fiterr]=fit_azi_anisotropy(azi,phv,weight);
    parastd=confint(para,.95);
    fit_azi.periods(ip)=periods(ip);
    fit_azi.c_iso(ip)=para.a;
    fit_azi.c_iso_95(ip)=parastd(2,1)-para.a;
    fit_azi.A2(ip)=para.d;
    fit_azi.A2_95(ip)=parastd(2,2)-para.d;
    fit_azi.phi2(ip)=para.e;
    fit_azi.phi2_95(ip)=parastd(2,3)-para.e;
end
 
fit_azi_bin = [];
weight = 0;
% Fit binned measurements
for ip = 1:length(periods)
    isgood = eventphv_ani(ip).isgood;
    phv = eventphv_ani(ip).phv(isgood);
    azi = eventphv_ani(ip).azi(isgood);
    azi(azi<0) = azi(azi<0)+360;
    w = eventphv_ani(ip).w(isgood);
    [~,Isort] = sort(azi);
    azi_srt = azi(Isort);
    phv_srt = phv(Isort);
    w_srt = w(Isort); %phv_std(Isort); % weight by std
    bins = [0:azi_bin_deg:360];
    phv_bin = 0; azi_bin = 0; phv_bin_err = 0; azi_bin_err = 0;
    for ibin = 1:length(bins)-1
        I_bin = azi_srt>=bins(ibin) & azi_srt<bins(ibin+1);
        if sum(I_bin)<min_nbin
            I_bin = false(size(I_bin));
        end
        wbin = w_srt(I_bin);
        wbinnorm = wbin/sum(wbin);
        % Weighted means and standard deviations
        fit_azi_bin.meas(ip).phv(ibin) = nansum(wbin .* phv_srt(I_bin)) / nansum(wbin);
        fit_azi_bin.meas(ip).phv_std(ibin) = sqrt( var(phv_srt(I_bin) , wbinnorm) );
        fit_azi_bin.meas(ip).azi(ibin) = (bins(ibin)+bins(ibin+1))/2;
        weight(ibin) = sum(I_bin);
    end
    
%     [para fiterr]=fit_azi_anisotropy(fit_azi_bin.meas(ip).azi, fit_azi_bin.meas(ip).phv, fit_azi_bin.meas(ip).phv_std.^(1/2) .* sqrt(weight));
%     [para fiterr]=fit_azi_anisotropy(fit_azi_bin.meas(ip).azi, fit_azi_bin.meas(ip).phv, sqrt(weight));
    [para fiterr]=fit_azi_anisotropy(fit_azi_bin.meas(ip).azi, fit_azi_bin.meas(ip).phv);
    parastd=confint(para,.95);
    fit_azi_bin.periods(ip)=periods(ip);
    fit_azi_bin.c_iso(ip)=para.a;
    fit_azi_bin.c_iso_95(ip)=parastd(2,1)-para.a;
    fit_azi_bin.A2(ip)=para.d;
    fit_azi_bin.A2_95(ip)=parastd(2,2)-para.d;
    fit_azi_bin.phi2(ip)=para.e;
    fit_azi_bin.phi2_95(ip)=parastd(2,3)-para.e;
end

fit_azi_bin_res = [];
weight = 0;
% Fit binned measurements
for ip = 1:length(periods)
    isgood = eventphv_ani(ip).isgood;
    phv = eventphv_ani(ip).phv(isgood);
    phv_iso = eventphv_ani(ip).phv_iso(isgood);
    dphv = (phv-phv_iso)./phv_iso;
    azi = eventphv_ani(ip).azi(isgood);
    azi(azi<0) = azi(azi<0)+360;
    w = eventphv_ani(ip).w(isgood);
    [~,Isort] = sort(azi);
    azi_srt = azi(Isort);
    dphv_srt = dphv(Isort);
    w_srt = w(Isort); %phv_std(Isort); % weight by std
    bins = [0:azi_bin_deg:360];
    phv_bin = 0; azi_bin = 0; phv_bin_err = 0; azi_bin_err = 0;
    for ibin = 1:length(bins)-1
        I_bin = azi_srt>=bins(ibin) & azi_srt<bins(ibin+1);
        if sum(I_bin)<min_nbin
            I_bin = false(size(I_bin));
        end
        wbin = w_srt(I_bin);
        wbinnorm = wbin/sum(wbin);
        % Weighted means and standard deviations
        fit_azi_bin_res.meas(ip).dphv(ibin) = nansum(wbin .* dphv_srt(I_bin)) / nansum(wbin);
        fit_azi_bin_res.meas(ip).dphv_std(ibin) = sqrt( var(dphv_srt(I_bin) , wbinnorm) );
        fit_azi_bin_res.meas(ip).azi(ibin) = (bins(ibin)+bins(ibin+1))/2;
        weight(ibin) = sum(I_bin);
    end
    
%     [para fiterr]=fit_azi_anisotropy(fit_azi_bin.meas(ip).azi, fit_azi_bin.meas(ip).phv, fit_azi_bin.meas(ip).phv_std.^(1/2) .* sqrt(weight));
%     [para fiterr]=fit_azi_anisotropy(fit_azi_bin.meas(ip).azi, fit_azi_bin.meas(ip).phv, sqrt(weight));
    [para fiterr]=fit_azi_anisotropy2theta_resid(fit_azi_bin_res.meas(ip).azi, fit_azi_bin_res.meas(ip).dphv);
    parastd=confint(para,.95);
    fit_azi_bin_res.periods(ip)=periods(ip);
    fit_azi_bin_res.A2(ip)=para.d;
    fit_azi_bin_res.A2_95(ip)=parastd(2,1)-para.d;
    fit_azi_bin_res.phi2(ip)=para.e;
    fit_azi_bin_res.phi2_95(ip)=parastd(2,2)-para.e;
end

matfilename = [eikonl_ani_output_path,'/eikonal_ani_',comp,'.mat'];
save(matfilename,'eventphv_ani','fit_azi','fit_azi_bin','fit_azi_bin_res');
disp(['Save the result to: ',matfilename])
 
 %%
	if isfigure
		N=3; M = floor(length(periods)/N) +1;
		figure(88)
		clf
		for ip = 1:length(periods)
			subplot(M,N,ip)
			ax = worldmap(lalim, lolim);
			set(ax, 'Visible', 'off')
			h1=surfacem(xi,yi,eventphv_ani(ip).GV);
			% set(h1,'facecolor','interp');
%			load pngcoastline
%			geoshow([S.Lat], [S.Lon], 'Color', 'black','linewidth',2)
			title(['Periods: ',num2str(periods(ip))],'fontsize',15)
			avgv = nanmean(eventphv_ani(ip).GV(:));
			if isnan(avgv)
				continue;
			end
			r = 0.1;
			caxis([avgv*(1-r) avgv*(1+r)])
			colorbar
			load seiscmap
			colormap(seiscmap)
            
%             % Plot ray paths
%             jj=0;
%             for ii = 1:length(eventphv_ani(ip).w)
%                 if ~eventphv_ani(ip).isgood(ii) %|| raytomo(ip).sta_dep(ixsp) > -4500
%                     continue
%                 end
%                 jj = jj + 1;
%                 lat1(jj,1) = eventphv_ani(ip).rays(ii,1);
%                 lon1(jj,1) = eventphv_ani(ip).rays(ii,2);
%                 lat2(jj,1) = eventphv_ani(ip).rays(ii,3);
%                 lon2(jj,1) = eventphv_ani(ip).rays(ii,4);
%             end
%             hold on;
%             h = plotm([lat1 lat2]',[lon1 lon2]','-k','linewidth',0.5); hold on;
		
        end
		drawnow;
	end
    
    %%
	figure(89);
	set(gcf,'position',[351   677   560   668]);
	clf
	clear avgv avgv_std aniso_str aniso_str_std aniso_azi aniso_azi_std
	for ip = 1:length(periods)
	    avgv(ip) = eventphv_ani(ip).GV_av;
	    avgv_std(ip) = 0; %nanmean(avgphv_aniso(ip).isophv_std(:));
	    aniso_str(ip) = eventphv_ani(ip).A2;
	    aniso_str_std(ip) = 0; %nanmean(avgphv_aniso(ip).aniso_strength_std(:));
	    aniso_azi(ip) = eventphv_ani(ip).phi2;
	    aniso_azi_std(ip) = 0; %nanmean(avgphv_aniso(ip).aniso_azi_std(:));
        
        c_iso(ip) = fit_azi.c_iso(ip);
        c_iso_95(ip) = fit_azi.c_iso_95(ip);
        A2(ip) = fit_azi.A2(ip);
        A2_95(ip) = fit_azi.A2_95(ip);
        phi2(ip) = fit_azi.phi2(ip);
        phi2_95(ip) = fit_azi.phi2_95(ip);
        
        c_iso_bin(ip) = fit_azi_bin.c_iso(ip);
        c_iso_95_bin(ip) = fit_azi_bin.c_iso_95(ip);
        A2_bin(ip) = fit_azi_bin.A2(ip);
        A2_95_bin(ip) = fit_azi_bin.A2_95(ip);
        phi2_bin(ip) = fit_azi_bin.phi2(ip);
        phi2_95_bin(ip) = fit_azi_bin.phi2_95(ip);
        
        A2_bin_res(ip) = fit_azi_bin_res.A2(ip);
        A2_95_bin_res(ip) = fit_azi_bin_res.A2_95(ip);
        phi2_bin_res(ip) = fit_azi_bin_res.phi2(ip);
        phi2_95_bin_res(ip) = fit_azi_bin_res.phi2_95(ip);
    end
    c = [3.9936 4.0041 4.0005 3.9999 3.9929 3.9832 3.9813 3.9841 3.9874 3.9996 4.0138 4.0519 4.0930 4.1677 4.2520];
    p2p = [0.0093 0.0119 0.0107 0.0081 0.0125 0.0138 0.0062 0.0027 0.0062 0.0090 0.0045 0.0060 0.0058 0.0051 0.0033];
    fastdir = [67.4502 76.4438 86.2857 82.7528 90.8716 93.6104 93.2020 106.5968 97.9209 102.4597 131.1436 123.1767 115.6864 101.7450 11.1088];
    pers = round(logspace(log10(20),log10(150),15));
    aniso_azi(aniso_azi<0) = aniso_azi(aniso_azi<0)+180;
	%plot native
	subplot(3,1,1); hold on;
    plot(pers,c,'-o','color',[0.8 0.8 0.8]);
	errorbar(periods,avgv,avgv_std*2,'-or');
    errorbar(periods,c_iso,c_iso_95,'-ob');
    errorbar(periods,c_iso_bin,c_iso_95_bin,'-ok');
	ylim([3.85 4.4]);
    xlim([20 150]);
	ylabel('c (km/s)');
	%plot native
	subplot(3,1,2); hold on;
    plot(pers,p2p*2*100,'-o','color',[0.8 0.8 0.8]);
	errorbar(periods,aniso_str*100*2,aniso_str_std*100*2,'-or');
    errorbar(periods,A2*100*2,A2_95*100,'-ob');
    errorbar(periods,A2_bin*100*2,A2_95_bin*100,'-ok');
    errorbar(periods,A2_bin_res*100*2,A2_95_bin_res*100,'-om');
	ylim([0 5]);
    xlim([20 150]);
	ylabel('2A');
	%plot native
	subplot(3,1,3);
	plot([periods(1),periods(end)],FSD*[1 1],'--k','linewidth',1.5); hold on;
	plot([periods(1),periods(end)],APM*[1 1],'--','color',[0.5 0.5 0.5],'linewidth',1.5);
	plot(pers,fastdir,'-o','color',[0.8 0.8 0.8]);
    errorbar(periods,aniso_azi,aniso_azi_std*2,'-or');
	errorbar(periods,aniso_azi+180,aniso_azi_std*2,'-or');
	errorbar(periods,aniso_azi-180,aniso_azi_std*2,'-or');
    errorbar(periods,phi2,phi2_95,'-ob');
    errorbar(periods,phi2+180,phi2_95,'-ob');
    errorbar(periods,phi2-180,phi2_95,'-ob');
    errorbar(periods,phi2_bin,phi2_95_bin,'-ok');
    errorbar(periods,phi2_bin+180,phi2_95_bin,'-ok');
    errorbar(periods,phi2_bin-180,phi2_95_bin,'-ok');
    errorbar(periods,phi2_bin_res,phi2_95_bin_res,'-om');
    errorbar(periods,phi2_bin_res+180,phi2_95_bin_res,'-om');
    errorbar(periods,phi2_bin_res-180,phi2_95_bin_res,'-om');
	ylim([50 180]);
    xlim([20 150]);
	ylabel('\phi');
	xlabel('Periods (s)');
    
    %%
    figure(90); clf;
    for ip = 1:length(periods)
        M=4;
        N=4;
        isgood = eventphv_ani(ip).isgood;
	    dt = eventphv_ani(ip).dt(isgood);
        azi = eventphv_ani(ip).azi(isgood);
        azi(azi<0) = azi(azi<0)+360;
        phv = eventphv_ani(ip).phv(isgood);
        avgv = eventphv_ani(ip).GV_av;
        phv_iso = eventphv_ani(ip).phv_iso(isgood);
        
        azi_bin = fit_azi_bin.meas(ip).azi;
        phv_bin = fit_azi_bin.meas(ip).phv;
        phv_iso_bin = fit_azi_bin.c_iso(ip);
        dv_std_bin = fit_azi_bin.meas(ip).phv_std./phv_iso_bin*100;
        dv_bin = (phv_bin-phv_iso_bin)./phv_iso_bin*100;
        azi_bin_res = fit_azi_bin_res.meas(ip).azi;
        dv_bin_res = fit_azi_bin_res.meas(ip).dphv*100;
        dv_std_bin_res = fit_azi_bin_res.meas(ip).dphv_std*100;
        
        dv = (phv-avgv)./avgv*100;
        dv2 = (phv-phv_iso)./phv_iso*100;
        dv2_not = (eventphv_ani(ip).ddist(~isgood)./eventphv_ani(ip).dt(~isgood)-eventphv_ani(ip).phv_iso(~isgood))./eventphv_ani(ip).phv_iso(~isgood)*100;
        dv3 = (phv-nanmean(phv_iso))./nanmean(phv_iso)*100;
        
        A2 = eventphv_ani(ip).A2;
        phi2 = eventphv_ani(ip).phi2;
        theta = [0:1:360];
        
        subplot(M,N,ip); hold on;
%         plot(azi,dv,'.'); hold on;
%         plot(eventphv_ani(ip).azi(~isgood),dv2_not,'.','color',[0.8 0.8 0.8]);
        plot(azi,dv2,'.','color',[0.8 0.8 0.8]); hold on;
        errorbar(azi_bin,dv_bin,dv_std_bin,'ok');
        errorbar(azi_bin_res,dv_bin_res,dv_std_bin_res,'om');
%         plot(azi,dv3,'.r');
        plot(theta,A2*cosd(2*(theta-phi2))*100,'-r');
        plot(theta,fit_azi.A2(ip)*cosd(2*(theta-fit_azi.phi2(ip)))*100,'-b');
        plot(theta,fit_azi_bin.A2(ip)*cosd(2*(theta-fit_azi_bin.phi2(ip)))*100,'-k');
        plot(theta,fit_azi_bin_res.A2(ip)*cosd(2*(theta-fit_azi_bin_res.phi2(ip)))*100,'-m');
        title([num2str(periods(ip)),' s']);
        ylim([-5 5]);
        xlim([0 360]);
    end
    
    figure(91); clf;
    for ip = 1:length(periods)
        isgood = eventphv_ani(ip).isgood;
	    dt = eventphv_ani(ip).dt(isgood);
        azi = eventphv_ani(ip).azi(isgood);
%         azi(azi<0) = azi(azi<0)+360;
        phv = eventphv_ani(ip).ddist(isgood)./dt;
        avgv = eventphv_ani(ip).GV_av;
        phv_iso = eventphv_ani(ip).phv_iso(isgood);
        dv = (phv-avgv)./avgv*100;
        dv2 = (phv-phv_iso)./phv_iso*100;
        
        A2 = eventphv_ani(ip).A2;
        phi2 = eventphv_ani(ip).phi2;
        subplot(M,N,ip);
        plot(azi,dv-A2*cosd(2*(azi-phi2))*100,'.'); hold on;
        plot(azi,dv2-A2*cosd(2*(azi-phi2))*100,'.g'); hold on;
%         ylim([-10 10]);
    end
