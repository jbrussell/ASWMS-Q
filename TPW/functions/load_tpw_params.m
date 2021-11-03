function [mat] = load_tpw_params(summfile,allfile)
%Load two-plane wave parameters (azimuth, init. amplitude, init. phase)
%

% Get number of events
cfp = fopen(allfile,'r');
stemp = fgetl(cfp);
eventnum = sscanf(stemp,'%d');
fclose(cfp);

% Read wave parameters and misfits
fid = fopen(summfile,'r');
% Skip header lines
for ii = 1:7
    fgetl(fid);
end
% Now loop through all events and get wave parameters
mat = [];
for ie = 1:eventnum
    % info
    stemp = fgetl(fid);
    dtemp = sscanf(stemp,' event            %d   %f  data std dev   %f rms phase misfit  s  rms amp mistfit      %f');
    mat.idnum(ie) = dtemp(1);
    mat.rmsdata(ie) = dtemp(2); % standard deviation of data
    mat.rmsphase(ie) = dtemp(3); % rms phase misfit (s)
    mat.rmsamp(ie) = dtemp(4); % rms amp misfit
    % wave 1
    stemp = fgetl(fid);
    dtemp = sscanf(stemp,'%f %f %f');
    mat.wvaz1(ie) = dtemp(1);
    mat.startamp1(ie) = dtemp(2);
    mat.stphase1(ie) = dtemp(3);
    % wave 2
    stemp = fgetl(fid);
    dtemp = sscanf(stemp,'%f %f %f');
    mat.wvaz2(ie) = dtemp(1);
    mat.startamp2(ie) = dtemp(2);
    mat.stphase2(ie) = dtemp(3);
end	
fclose(fid);

end

