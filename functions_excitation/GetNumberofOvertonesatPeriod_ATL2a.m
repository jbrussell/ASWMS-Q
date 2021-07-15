function [ maxN ] = GetNumberofOvertonesatPeriod_ATL2a(PeriodQ )
% how many overtones exist at this period?

info = load('OvertoneNumber_MaxPeriod.txt');
nvals=info(:,1);
periodlist=info(:,2);


Vq = interp1(periodlist,nvals,PeriodQ);
maxN=floor(Vq);
end

