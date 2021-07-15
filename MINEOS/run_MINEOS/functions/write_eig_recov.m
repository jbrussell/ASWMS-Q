% 10/7/15 -- Josh Russell
%
% Function to write the input file for eig_recover
%

function write_eig_recov(FNUM,ll)
parameter_FRECHET;
TYPE = param.TYPE;

if strcmp(TYPE,'S') == 1
    TYPEID = param.STYPEID;
    RUNFILE_eigrecov = 'run_eigrecov.s';
elseif strcmp(TYPE,'T') == 1
    TYPEID = param.TTYPEID;
    RUNFILE_eigrecov = 'run_eigrecov.t';
end

EIG = [CARDTABLE,param.CARDID,'.',TYPEID,'_',num2str(FNUM),'.eig'];

% Write eig_recover driver
fid=fopen(RUNFILE_eigrecov,'w');
fprintf(fid,'%s\n',EIG);
fprintf(fid,'%d\n',ll);
fprintf(fid,'\n');
fclose(fid);
