%% Write run file for frechet.f
% NJA, 2014
% make TYPEID as a parameter in parameter_FRECHET
% pylin.patty 2015/01
%
function write_frechet(TYPE,CARDID,NDISC,ZDISC)

parameter_FRECHET;
CARDPATH  = param.CARDPATH;
TABLEPATH = param.TABLEPATH;
% FRECHETPATH = param.frechetpath;
FRECHETPATH = [param.frechet,CARDID,'/'];
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
    EIG = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.eig'];

    
elseif strcmp(TYPE,'S') == 1
    disp('Spheroidal!');
    
    RUNFILE = 'run_frechet.s';
    TYPEID = param.STYPEID;
    
else
    disp('No TYPE recognized!');
    
end

QMOD = [CARDPATH,CARDID,'.qmod'];
BRANCH = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.table_hdr.branch'];
FRECH = [FRECHETPATH,CARDID,'.',TYPEID,'.frech'];

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

if strcmp(TYPE,'T') == 1
    fprintf(fid,'%s\n',EIG);
    fprintf(fid,'%i\n',NDISC);
    
    
    for idisc = 1:NDISC
        fprintf(fid,'%i\n',ZDISC(idisc));
    end
    fprintf(fid,'\n');
    fclose(fid);
    
elseif strcmp(TYPE,'S') == 1
    %FREQ = '0_50';
    %EIG1 = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.',FREQ,'.eig'];    
    EIG1 = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.eig'];
    
    fprintf(fid,'%s\n',EIG1);
    fprintf(fid,'%i\n',NDISC);
    
    for idisc = 1:NDISC
        fprintf(fid,'%i\n',ZDISC(idisc));
    end

    
    fprintf(fid,'\n');
    fclose(fid);
    
end
