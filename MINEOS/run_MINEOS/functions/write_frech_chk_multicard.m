% 10/7/15 -- Josh Russell
%
% Creates the input file for frechet when there are mutiple eigenfunction
% files.
%
% 11/2/15 JBR -- This version allows for custom CARDID input
%
function write_frech_chk_multicard(NDISC,CARDID)


parameter_FRECHET;
CARDPATH  = param.CARDPATH;
TABLEPATH = param.TABLEPATH;
% FRECHETPATH = param.frechetpath;
FRECHETPATH = [param.frechet,CARDID,'/'];
TYPE = param.TYPE;
% Set number of discontinuities to add
if NDISC > 0
    
    if length(ZDISC) > NDISC
        disp('Mismatch in discontinuity depths!')
    end
end

if strcmp(TYPE,'T') == 1
    disp('Toroidal!');
    
    RUNFILE = 'run_frechet.t';
    TYPEID = param.TTYPEID;

    
elseif strcmp(TYPE,'S') == 1
    disp('Spheroidal!');
    
    RUNFILE = 'run_frechet.s';
    TYPEID = param.STYPEID;
    
else
    disp('No TYPE recognized!');
    
end

EIG0 = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'_0.eig_fix'];
QMOD = [CARDPATH,CARDID,'.qmod'];
BRANCH = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.table_hdr.branch'];
FRECH = [FRECHETPATH,CARDID,'.',TYPEID,'.frech'];

com = ['ls ',TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'_*.eig_fix | cat'];
[status eig_fils] = system(com);
EIG = strsplit(eig_fils,'\n');

% Check to see if file exists ... program will not overwrite it if it does
if exist(FRECH,'file') == 2
    disp([FRECH,' File exists! Removing it now'])
    com = ['rm -f ',FRECH];
    [status,log] = system(com);
end


fid = fopen(RUNFILE,'w');
fprintf(fid,'%s\n',QMOD);
fprintf(fid,'%s\n',BRANCH);
fprintf(fid,'%s\n',FRECH);
fprintf(fid,'%s\n',EIG0);
fprintf(fid,'%i\n',NDISC);
for i = 2:size(EIG,2)-1 % skip *_0.eig_fix and start at *_1.eig_fix
    disp(['Using eig file ',num2str(i-1)])
    fprintf(fid,'%s\n',EIG{i});
end
fprintf(fid,'\n');
fclose(fid);
    

    