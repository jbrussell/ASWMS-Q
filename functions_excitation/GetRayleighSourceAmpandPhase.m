function [ SourceAmp,SourcePhase,Complex_Rad_Pattern,Term1,Term2,Term3 ] = GetRayleighSourceAmpandPhase(AZI,...
    Source_Depth_m,period,Radius_List_m,ULIST,UderivLIST,VLIST,...
    VderivLIST,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp,Phvel,wavegroup_index )
% Convert the complex radiation pattern into meaningful phase and
% amplitude, as a function of azimuth...

% [ Complex_Rad_Pattern,Term1,Term2,Term3 ] = Get_Rayleigh_Rad_Pattern(AZI,Source_Depth_m,period,...
% Radius_List_m,ULIST,UderivLIST,VLIST,VderivLIST,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp,Phvel );

% wavegroup_index s = 1 minor arc or 2 major arc 
% 
 [ Complex_Rad_Pattern,Term1,Term2,Term3 ] = Get_Rayleigh_Rad_Pattern_DahlenTromp(AZI,Source_Depth_m,period,...
 Radius_List_m,ULIST,UderivLIST,VLIST,VderivLIST,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp,Phvel,wavegroup_index );

SourceAmp = abs(Complex_Rad_Pattern);
Real_Part = real(Complex_Rad_Pattern);
Imag_Part = imag(Complex_Rad_Pattern);
SourcePhasedecoy = atan2(Imag_Part,Real_Part);
%SourcePhase = atan(Imag_Part./Real_Part);
SourcePhase = zeros(size(SourcePhasedecoy));
% 
 SourcePhase(SourcePhasedecoy > 0) = SourcePhasedecoy(SourcePhasedecoy > 0)-pi;
SourcePhase(SourcePhasedecoy < 0) = SourcePhasedecoy(SourcePhasedecoy < 0)+pi;

end
