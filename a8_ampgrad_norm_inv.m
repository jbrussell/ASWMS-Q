% Read in the eventcs structures and apply eikonal tomography on each event.
% Written by Ge Jin, jinwar@gmail.com
% 2013.1.16
%
% jbrussell - modify original eikonal tomography which inverts dt
% measurements for dt/dx = grad(t) to instead invert dA/A measurements for
% (dA/dx)/A = grad(A)/A
%
clear
% setup parameters
setup_parameters

% debug setting
isfigure = 0;
isdisp = 0;
is_overwrite = 1;

is_receiver_terms = parameters.is_receiver_terms; % Correct amplitudes using receiver terms calculated from a8_0_receiver_terms?
is_offgc_smoothing = parameters.is_offgc_smoothing; % allows off-great-circle 1st derivative smoothing. Requires an initial run of a6_a0_eikonal_eq_GetPropAzi.m to get propagation azimuth

workingdir = parameters.workingdir;
receiverterms_path = [workingdir];
% input path
eventcs_path = [workingdir,'CSmeasure/'];
eikonl_propazi_output_path = [workingdir,'eikonal_propazi/'];
% output path
ampgrad_norm_output_path = [workingdir,'ampgrad_norm/'];

if ~exist(ampgrad_norm_output_path)
	mkdir(ampgrad_norm_output_path);
end

comp = parameters.component;
lalim=parameters.lalim;
lolim=parameters.lolim;
gridsize=parameters.gridsize;
periods = parameters.periods;
raydensetol=parameters.raydensetol;
smweight_array = parameters.smweight_array;
flweight_array = parameters.flweight_array; % JBR
Tdumpweight0 = parameters.Tdumpweight;
Rdumpweight0 = parameters.Rdumpweight;
fiterrtol = parameters.fiterrtol;
dterrtol = parameters.dterrtol;
isRsmooth = parameters.isRsmooth;
inverse_err_tol = parameters.inverse_err_tol;
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
    % dy = diff(ynode)*cosd(mean(xnode));  % correct smoothing for latitude
    % dy1 = dy(j(:)-1);
    % dy2 = dy(j(:));
    dy1 = km2deg(distance(xnode(i(:)),ynode(j(:)),xnode(i(:)),ynode(j(:)-1),referenceEllipsoid('GRS80'))/1000);
    dy2 = km2deg(distance(xnode(i(:)),ynode(j(:)),xnode(i(:)),ynode(j(:)+1),referenceEllipsoid('GRS80'))/1000);

    Areg = sparse(repmat(ind,1,3),[ind-1,ind,ind+1], ...
                    [-2./(dy1.*(dy1+dy2)), 2./(dy1.*dy2), -2./(dy2.*(dy1+dy2))],Nx*Ny,Nx*Ny);

    % latitude smoothing
    [i,j] = ndgrid(2:(Nx-1),1:Ny);
    ind = j(:) + Ny*(i(:)-1);
    % dx = diff(xnode);
    % dx1 = dx(i(:)-1);
    % dx2 = dx(i(:));
    dx1 = km2deg(distance(xnode(i(:)),ynode(j(:)),xnode(i(:)-1),ynode(j(:)),referenceEllipsoid('GRS80'))/1000);
    dx2 = km2deg(distance(xnode(i(:)),ynode(j(:)),xnode(i(:)+1),ynode(j(:)),referenceEllipsoid('GRS80'))/1000);

    Areg = [Areg;sparse(repmat(ind,1,3),[ind-Ny,ind,ind+Ny], ...
            [-2./(dx1.*(dx1+dx2)), 2./(dx1.*dx2), -2./(dx2.*(dx1+dx2))],Nx*Ny,Nx*Ny)];

    F=sparse(Nx*Ny*2*2,Nx*Ny*2);
    for n=1:size(Areg,1)
        ind=find(Areg(n,:)~=0);
        F(2*n-1,2*ind-1)=Areg(n,ind);
        F(2*n,2*ind)=Areg(n,ind);
    end
toc

% JBR - define first derivative "flatness" kernel
F2 = flat_kernel_build(xnode, ynode, Nx*Ny);

% read in bad station list, if existed
if exist('badampsta.lst')
	badstnms = textread('badampsta.lst','%s');
	disp('Found Bad stations:')
	disp(badstnms)
end

if is_receiver_terms==1
    load([receiverterms_path,'receiver_terms_',parameters.component,'.mat']);
end

csmatfiles = dir([eventcs_path,'/*cs_',comp,'.mat']);
for ie = 1:length(csmatfiles)
%for ie = 30
	clear ampgrad_norm 
	% read in data and set up useful variables
	temp = load([eventcs_path,csmatfiles(ie).name]);
	eventcs =  temp.eventcs;
	disp(eventcs.id)
	evla = eventcs.evla;
	evlo = eventcs.evlo;

	matfilename = [ampgrad_norm_output_path,'/',eventcs.id,'_eikonal_',comp,'.mat'];
	if exist(matfilename,'file') && ~is_overwrite
		disp(['Exist ',matfilename,', skip']);
		continue;
	end


	if exist('badstnms','var')
		badstaids = find(ismember(eventcs.stnms,badstnms));
	else
		badstaids = [];
	end

	% Load previous eikonal mat to get propagation azimuth
	if is_offgc_smoothing==1
	    eikonal_in = [eikonl_propazi_output_path,'/',eventcs.id,'_eikonal_',comp,'.mat'];
	    if ~exist(eikonal_in,'file')
	        error('No propagation azimuth found. Need to first run a6_a0_eikonal_eq_GetPropAzi.m');
	    end
	    aziprop = load(eikonal_in);
	end

	% Build the ray locations
	clear rays 
	for ics = 1:length(eventcs.CS)
		rays(ics,1) = eventcs.stlas(eventcs.CS(ics).sta1);
		rays(ics,2) = eventcs.stlos(eventcs.CS(ics).sta1);
		rays(ics,3) = eventcs.stlas(eventcs.CS(ics).sta2);
		rays(ics,4) = eventcs.stlos(eventcs.CS(ics).sta2);
	end

	% Build the kernel
	disp('Buildling up ray path kernel')
	tic
		mat=kernel_build(rays,xnode,ynode);
	toc

	% Loop through the periods
	for ip = 1:length(periods)
		
		% Build the rotation matrix
	    if is_offgc_smoothing==1
	        phase_lat = -aziprop.eventphv(ip).GVx; % phase slowness in x-direction
	        phase_lon = -aziprop.eventphv(ip).GVy; % phase slowness in y-direction
	        razi = 90 - atan2d(phase_lat,phase_lon);
	        azimat_ev = azimuth(xi+gridsize/2,yi+gridsize/2,evla,evlo,referenceEllipsoid('GRS80'))+180;
	        razi(isnan(razi)) = azimat_ev(isnan(razi));
	    else
	        razi = azimuth(xi+gridsize/2,yi+gridsize/2,evla,evlo,referenceEllipsoid('GRS80'))+180;
	    end
	    R = sparse(2*Nx*Ny,2*Nx*Ny);
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
		
		% build dumping matrix for ST
	    dumpmatT = R(2:2:2*Nx*Ny,:);

	    % build dumping matrix for SR
	    dumpmatR = R(1:2:2*Nx*Ny-1,:);
		
		smweight0 = smweight_array(ip);
		flweight0 = flweight_array(ip); % JBR
		dt = zeros(length(eventcs.CS),1);
	    
	    % Load amplitude data
	    amps = zeros(1,length(eventcs.stlas));
	    dist = zeros(1,length(eventcs.stlas));
		for ista = 1:length(eventcs.autocor)
	        dist(ista) = vdist(eventcs.evla,eventcs.evlo,eventcs.stlas(ista),eventcs.stlos(ista))/1000;
			if eventcs.autocor(ista).exitflag(ip)>0
				amps(ista) = eventcs.autocor(ista).amp(ip);
			else
				amps(ista) = NaN;
			end
		end
		% change from power spectrum to amplitude
		amps = amps.^.5;
	    
	    % Correct amplitude for local receiver effects
	    if is_receiver_terms==1
	        Amp_rec = receiver(ip).Amp_rec;
	        for ista = 1:length(eventcs.stnms)
	            Istation = find(strcmp(eventcs.stnms(ista),receiver(ip).stas));
	            if isempty(Istation)
	                disp(['No station term for ',eventcs.stnms(ista)]);
	                continue
	            end
	            amps(ista) = amps(ista) ./ Amp_rec(Istation);
	        end
	    end
	    
	    
	    amplitudes(ip).amps = amps;
	    
	    
		w = zeros(length(eventcs.CS),1);
		for ics = 1:length(eventcs.CS)
	        Ista1 = eventcs.CS(ics).sta1;
	        Ista2 = eventcs.CS(ics).sta2;
			if eventcs.CS(ics).isgood(ip) > 0 && ~isnan(amps(Ista1)*amps(Ista2))
	            
	%                 dt(ics) = amps(Ista1)-amps(Ista2);
	            amp_mean = 0.5*(amps(Ista1)+amps(Ista2));
	            dt(ics) = (amps(Ista1)-amps(Ista2)) ./ amp_mean;
	            
	% 				dt(ics) = eventcs.CS(ics).dtp(ip);
				w(ics) = 1;
			else
				dt(ics) = 0;
				w(ics) = 0;
			end
			if sum(ismember([eventcs.CS(ics).sta1 eventcs.CS(ics).sta2],badstaids)) > 0
				w(ics) = 0;
			end
		end
		W = sparse(length(w),length(w));
		for id = 1:length(w)
	        if w(id) > 0
	            W(id,id) = w(id);
	        end
	    end

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

		% Set up matrix on both side
		if isRsmooth
	        A=[W*mat;smweight*F*R;flweight*F2*R;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
	    else
	        A=[W*mat;smweight*F;flweight*F2;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
	    end

		avgv = eventcs.avgphv(ip);
	    rhs=[W*dt;zeros(size(F,1),1);zeros(size(F2,1),1);zeros(size(dumpmatT,1),1);dumpweightR*ones(size(dumpmatR,1),1)./avgv];
	    
		% Least square inversion
		if isempty(W(W~=0)) || ~isempty(W(isnan(W))) || ~isempty(W(isinf(W)))
            % Skip if no good data or if W is nan
            disp('No good data or NaNs in W matrix, skipping...');
            phaseg = nan(size(A,2),1);
            A = eye(size(A));
        else
            phaseg=(A'*A)\(A'*rhs);
        end
	        
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
	        err = mat*phaseg - dt;
			err(diag(W)==0) = 0;
	        stderr=std(err(err~=0));
	        if stderr > dterrtol
	            stderr = dterrtol;
	        end
	        for i=1:length(err)
	            if abs(err(i)) > inverse_err_tol*stderr  || abs(err(i))==0
	                W(i,i)=0;
				else
					W(i,i)=1./stderr;
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
	        
			if isRsmooth
	            A=[W*mat;smweight*F*R;flweight*F2*R;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
	        else
	            A=[W*mat;smweight*F;flweight*F2;dumpweightT*dumpmatT;dumpweightR*dumpmatR];
	        end
	        rhs=[W*dt;zeros(size(F,1),1);zeros(size(F2,1),1);zeros(size(dumpmatT,1),1);dumpweightR*ones(size(dumpmatR,1),1)./avgv];
			if isempty(W(W~=0)) || ~isempty(W(isnan(W))) || ~isempty(W(isinf(W)))
	            % Skip if no good data or if W is nan
	            disp('No good data or NaNs in W matrix, skipping...');
	            phaseg = nan(size(A,2),1);
	            A = eye(size(A));
	        else
	            phaseg=(A'*A)\(A'*rhs);
	        end
	    end	

		% Estimate differential amplitude residuals
	    dA_res = dt - mat*phaseg;
		
		% Calculate model resolution and chi2
	    Ginv = (A'*A)\mat';
	    R = Ginv * mat; % model resolution
	    D = mat * Ginv; % data resolution
	    % Effective degrees of freedom
	    v = length(dt) - trace(D);
	%         v = trace(D);
	    % normalized chi2 uncertainties
	    res = (dt-mat*phaseg);
	    res(diag(W)==0) = nan;
	    rms_res = sqrt(nanmean(res.^2));
	    dt_std = rms_res;
	    chi2 = nansum(res.^2./dt_std.^2)/v;

	    % Calculate model uncertainties
	    dAmp_std = diag(inv(A'*A)).^(1/2);
	    for i=1:Nx
			for j=1:Ny
				n=Ny*(i-1)+j;
				dAmpx_err(i,j)= dAmp_std(2*n-1);
				dAmpy_err(i,j)= dAmp_std(2*n);
			end
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

		% Change phaseg into phase velocity
		for i=1:Nx
			for j=1:Ny
				n=Ny*(i-1)+j;
				dAmpx(i,j)= phaseg(2*n-1);
				dAmpy(i,j)= phaseg(2*n);
			end
		end
		dAmp=(dAmpx.^2+dAmpy.^2).^.5;
		% Get rid of uncertain area
		dAmpx_err(isnan(dAmp)) = nan;
	    dAmpy_err(isnan(dAmp)) = nan;
	    
	    % Propagate errors
	    dAmp_err = (((dAmpx.*dAmpx_err).^2 + (dAmpy.*dAmpy_err).^2)).^0.5 ./ (dAmpx.^2+dAmpy.^2).^0.5; % propagate errors dAmp magnitude

		% save the result in the structure
		ampgrad_norm(ip).rays = rays;
		ampgrad_norm(ip).w = diag(W);
		ampgrad_norm(ip).goodnum = length(find(w>0));
		ampgrad_norm(ip).badnum = length(find(w==0));
		ampgrad_norm(ip).dt = dt;
		ampgrad_norm(ip).dA_res = dA_res; % data residuals
		ampgrad_norm(ip).chi2 = chi2; % chi2 misfit
		ampgrad_norm(ip).dAmp_A = dAmp;
		ampgrad_norm(ip).dAmpx_A = dAmpx;
		ampgrad_norm(ip).dAmpy_A = dAmpy;
		ampgrad_norm(ip).dAmp_A_err = dAmp_err;
		ampgrad_norm(ip).dAmpx_A_err = dAmpx_err;
		ampgrad_norm(ip).dAmpy_A_err = dAmpy_err;
		ampgrad_norm(ip).raydense = raydense;
		ampgrad_norm(ip).lalim = lalim;
		ampgrad_norm(ip).lolim = lolim;
		ampgrad_norm(ip).gridsize = gridsize;
		ampgrad_norm(ip).id = eventcs.id;
		ampgrad_norm(ip).evla = eventcs.evla;
		ampgrad_norm(ip).evlo = eventcs.evlo;
		ampgrad_norm(ip).evdp = eventcs.evdp;
		ampgrad_norm(ip).period = periods(ip);
		ampgrad_norm(ip).amps = amplitudes(ip).amps;
		ampgrad_norm(ip).stlas = eventcs.stlas;
		ampgrad_norm(ip).stlos = eventcs.stlos;
		ampgrad_norm(ip).stnms = eventcs.stnms;
		ampgrad_norm(ip).isgood = ampgrad_norm(ip).w>0;
		ampgrad_norm(ip).Mw = eventcs.Mw;
		disp(['Period:',num2str(periods(ip)),', Goodnum:',num2str(ampgrad_norm(ip).goodnum),...
				'Badnum:',num2str(ampgrad_norm(ip).badnum)]);
	end % end of periods loop
	if isfigure
		N=3; M = floor(length(periods)/N) +1;
		figure(88)
		clf
	    sgtitle('Amplitude gradient','fontsize',18,'fontweight','bold');
	    lims = {};
		for ip = 1:length(periods)
			subplot(M,N,ip)
			ax = worldmap(lalim, lolim);
			set(ax, 'Visible', 'off')
			h1=surfacem(xi,yi,ampgrad_norm(ip).dAmp_A);
			% set(h1,'facecolor','interp');
	%			load pngcoastline
	%			geoshow([S.Lat], [S.Lon], 'Color', 'black','linewidth',2)
			title(['Periods: ',num2str(periods(ip))],'fontsize',15)
			avgv = nanmean(ampgrad_norm(ip).dAmp_A(:));
			if isnan(avgv)
				continue;
			end
			r = 0.8;
			caxis([avgv*(1-r) avgv*(1+r)])
			cb = colorbar;
	        lims{ip} = cb.Limits;
			load seiscmap
			colormap(seiscmap)
		end
		drawnow;
	    
	%         figure(87)
	% 		clf
	%         sgtitle('Amplitude gradient (surface calculation)','fontsize',18,'fontweight','bold');
	% 		for ip = 1:length(periods)
	%             attenuation_path = ['NoMelt_M5.5_detrend_Zcorr_100km_snr3_600km_SAVE2_nostationterms/attenuation/',eventcs.id,'_attenuation_BHZ.mat'];
	%             temp = load(attenuation_path);
	%             attenuation = temp.attenuation;
	%             amp_grad = attenuation(ip).amp_grad';
	%             inan = isnan(ampgrad_norm(ip).dAmp_A);
	%             amp_grad(inan) = nan;
	% 			subplot(M,N,ip)
	% 			ax = worldmap(lalim, lolim);
	% 			set(ax, 'Visible', 'off')
	% 			h1=surfacem(xi,yi,amp_grad);
	% 			% set(h1,'facecolor','interp');
	% %			load pngcoastline
	% %			geoshow([S.Lat], [S.Lon], 'Color', 'black','linewidth',2)
	% 			title(['Periods: ',num2str(periods(ip))],'fontsize',15)
	% 			if isnan(avgv)
	% 				continue;
	% 			end
	% 			colorbar
	%             if ~isempty(lims{ip})
	%                 caxis(lims{ip});
	%             end
	% 			load seiscmap
	% 			colormap(seiscmap)
	% 		end
	% 		drawnow;
	end
	matfilename = [ampgrad_norm_output_path,'/',eventcs.id,'_ampgrad_norm_',comp,'.mat'];
	save(matfilename,'ampgrad_norm');
	disp(['Save the result to: ',matfilename])
end % end of loop ie

%%% Plot residuals
eventfiles = dir([ampgrad_norm_output_path,'/*_ampgrad_norm_',parameters.component,'.mat']);
clear residuals
for ie = 1:length(eventfiles)
    temp = load(fullfile(ampgrad_norm_output_path,eventfiles(ie).name));
    ampgrad_norm = temp.ampgrad_norm;
    for ip = 1:length(ampgrad_norm)
        if ie == 1
            residuals(ip).rms_dA_res = [];
            residuals(ip).mean_dA_res = [];
        end
        isgood = ampgrad_norm(ip).isgood;
        dA_res = ampgrad_norm(ip).dA_res(isgood);
        residuals(ip).rms_dA_res = [residuals(ip).rms_dA_res(:); rms(dA_res(:))];
        residuals(ip).mean_dA_res = [residuals(ip).mean_dA_res(:); mean(dA_res(:))];
    end
end

%%
figure(86); clf; set(gcf,'color','w','position',[1035         155         560         781]);
for ip = 1:length(periods)
    subplot(2,1,1);
    plot(periods(ip),residuals(ip).mean_dA_res,'o','color',[0.7 0.7 0.7]); hold on;
    plot(periods(ip),nanmean(residuals(ip).mean_dA_res),'rs','linewidth',2,'markersize',10);
    ylabel('mean (dA/A_{obs}-dA/A_{pre})')
    set(gca,'linewidth',1.5,'fontsize',15);
    subplot(2,1,2);
    plot(periods(ip),residuals(ip).rms_dA_res,'o','color',[0.7 0.7 0.7]); hold on;
    plot(periods(ip),nanmean(residuals(ip).rms_dA_res),'rs','linewidth',2,'markersize',10);
    xlabel('Period (s)');
    ylabel('RMS (dA/A_{obs}-dA/A_{pre})')
    set(gca,'linewidth',1.5,'fontsize',15);
end
