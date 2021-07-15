% 10/7/15 -- Josh Russell
%
% Function to build the Mode.in files while looping through the mineos
% check
%
% INPUT: l_start - Starting l
%        LOOP - the loop number for the mineos_nohang run
%
% OUTPUT: Prints mode file to .../MODE_tables/MODE.in
%

function write_mode_in(l_start,LOOP)

parameter_FRECHET
TYPE = param.TYPE;
if strcmp(TYPE,'S') == 1
    if isnan(LOOP)
        MODEOUT = [param.MODEPATH,param.SMODEIN];
    else
        MODEOUT = [param.MODEPATH,param.SMODEIN,'_',num2str(LOOP)];
    end
    MODETYPE = 3;
elseif strcmp(TYPE,'T') == 1
    if isnan(LOOP)
        MODEOUT = [param.MODEPATH,param.TMODEIN];
    else
        MODEOUT = [param.MODEPATH,param.TMODEIN,'_',num2str(LOOP)];
    end
    MODETYPE = 2;
end
fid = fopen(MODEOUT,'w');
fprintf(fid,'1.d-12  1.d-12  1.d-12 .126\n');
fprintf(fid,'%d\n',MODETYPE);
fprintf(fid,'%d %d %.2f %.2f %d\n',[l_start maxL minF+0.05 maxF+0.05 N_modes]);
fprintf(fid,'0\n');
fclose(fid);