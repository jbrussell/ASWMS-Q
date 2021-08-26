
function [eventcs] = GetExcitationRatios(eventcs,Periodlist,CARD,Path2Phv,Path2Eig,MaxN,MinorOrMajor,comp)

if contains(comp,'Z') || contains(comp,'R') || contains(comp,'BDH') % Rayleigh
    eventcs = Rayleigh_GetExcitationRatios(eventcs,Periodlist,CARD,Path2Phv,Path2Eig,MaxN,MinorOrMajor);
elseif contains(comp,'T') % Love
    eventcs = Love_GetExcitationRatios(eventcs,Periodlist,CARD,Path2Phv,Path2Eig,MaxN,MinorOrMajor);
end
