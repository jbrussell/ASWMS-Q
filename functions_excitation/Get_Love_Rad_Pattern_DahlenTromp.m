function [ Complex_Rad_Pattern,Term1,Term2 ] = Get_Love_Rad_Pattern_DahlenTromp(AZI,Source_Depth_m,...
    period,Radius_List_m,WLIST,WderivLIST,...
    Mrr,Mtt,Mpp,Mrt,Mrp,Mtp,Love_Phvel,wavegroup_index )
% Get the (complex) radiation pattern for Love waves by taking in the inputs
% of W eignfunctions and their derivatives wrt depth.
% Also takes in a moment tensor and components
% Uses equation 25 in Wang and Dahlen
% Assumes the eigenfunction is very well sampled in radius space
% Assumes the depth you're interested in is in METERS

RE_m = 6371000;
omega = 2*pi/period;
k = omega/Love_Phvel;
AZI=Map_ClockwizeAzi2CounterClockwisefromSAzi( AZI );
s=wavegroup_index;

% Get index corresponding to your source depth of interest
source_rad = RE_m-Source_Depth_m;
difference_list = abs(source_rad-Radius_List_m);
[meaninglessval,mindx] = min(difference_list);
% normalize
source_rad = Radius_List_m(mindx)/RE_m;
%source_rad = Radius_List_m(mindx);

% Get the value of the derivatives of W, and W itself, here
Wdot = WderivLIST(mindx);
W = WLIST(mindx);

% Dahlen and Tromp Version
Term1 = -1.*(1./source_rad).*k.*W.*(-1*Mtp.*cosd(2.*AZI) +...
    0.5.*(Mtt-Mpp).*sind(2.*AZI));
Term2 = ((-1)^s).*(Wdot-W/source_rad).*(Mrt.*sind(AZI) - Mrp.*cosd(AZI));

Complex_Rad_Pattern = omega*(Term1.*exp(1i*pi/4)+Term2.*exp(-1*1i*pi/4));
%Complex_Rad_Pattern = omega*(Term3);

end

