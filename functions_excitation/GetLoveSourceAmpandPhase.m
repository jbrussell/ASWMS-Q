function [ SourceAmp,SourcePhase,Complex_Rad_Pattern,Term1,Term2 ] = GetLoveSourceAmpandPhase(AZI,...
    Source_Depth_m,period,Radius_List_m,WLIST,WderivLIST,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp,Love_Phvel,wavegroup_index  )
% Convert the complex radiation pattern into meaningful phase and
% amplitude, as a function of azimuth...

[ Complex_Rad_Pattern,Term1,Term2 ] = Get_Love_Rad_Pattern_DahlenTromp(AZI,...
Source_Depth_m,period,Radius_List_m,WLIST,WderivLIST,...
Mrr,Mtt,Mpp,Mrt,Mrp,Mtp,Love_Phvel,wavegroup_index  );

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
