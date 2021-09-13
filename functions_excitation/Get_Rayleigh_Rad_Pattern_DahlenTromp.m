function [ Complex_Rad_Pattern,Term1,Term2,Term3 ] = Get_Rayleigh_Rad_Pattern_DahlenTromp(AZI,...
    Source_Depth_m,period,Radius_List_m,ULIST,UderivLIST,VLIST,VderivLIST,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp,Phvel,wavegroup_index)

% Get the (complex) radiation pattern by taking in the inputs
% of U, V eignfunctions and their derivatives wrt depth.
% Also takes in a moment tensor and components
% Uses equation 11.34 in Dahlen and Tromp
% Assumes the eigenfunction is very well sampled in radius space
% Assumes the depth you're interested in is in METERS

RE_m = 6371000;
omega = 2*pi/period;
k = omega/Phvel;
AZI=Map_ClockwizeAzi2CounterClockwisefromSAzi( AZI );

% Get index corresponding to your source depth of interest
source_rad = RE_m-Source_Depth_m;
difference_list = abs(source_rad-Radius_List_m);
[meaninglessval,mindx] = min(difference_list);
% normalize
source_rad = Radius_List_m(mindx)/RE_m;
%source_rad = Radius_List_m(mindx);

% Get the value of the derivatives of U and V AND U AND V here
Udot = UderivLIST(mindx);
Vdot = VderivLIST(mindx);
U = ULIST(mindx);
V = VLIST(mindx);

s = wavegroup_index;

% Dahlen and tromp Version
Term1 = Udot*Mrr + (Mtt+Mpp).*(1/source_rad).*(U-0.5*k*V);
Term2 = -1*(k*V./source_rad).*(Mtp.*sind(2.*AZI ) +  0.5.*(Mtt-Mpp).* cosd(2.*AZI ));
Term3 = ((-1)^s)*(Vdot+(1./source_rad).*(k*U-V)).*(Mrp.*sind(AZI ) + Mrt.*cosd(AZI ));

% Complex_Rad_Pattern = omega*(Term1.*exp(1i*pi/4)+Term2.*exp(-1*1i*pi/4)+Term3.*exp(1i*pi/4));
Complex_Rad_Pattern = omega*(Term1.*exp(1i*pi/4)+Term3.*exp(-1*1i*pi/4)+Term2.*exp(1i*pi/4));
%Complex_Rad_Pattern = omega*(Term3);

% JBR scale to displacement
idx_surf = max(find(ULIST~=0));
U_surf = ULIST(idx_surf);
Complex_Rad_Pattern = Complex_Rad_Pattern .* (U_surf);

end
