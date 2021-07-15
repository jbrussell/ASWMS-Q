%% Write idagrn driver
% JBR 7/18
%
function write_idagrn(TYPE,CARD,EVT,STA,LENGTH_HR,DT,COMP)


parameter_FRECHET;

CARDTABLE = CARDTABLE;
CARDID = param.CARDID;

CARDPATH = param.CARDPATH;
FRECHETPATH = param.frechetpath;
TABLEPATH = param.TABLEPATH;

if strcmp(TYPE,'T') == 1
    disp('Toroidal!');
    
    RUNFILE = 'run_idagrn.t';
    TYPEID = param.TTYPEID;
    
elseif strcmp(TYPE,'S') == 1
    disp('Spheroidal!');
    
    RUNFILE = 'run_idagrn.s';
    TYPEID = param.STYPEID;
    
else
    disp('No TYPE recognized!');
    
end
TABLE = [CARDTABLE,CARD,'.',TYPEID,'.table'];

%% Get event information
com = sprintf('awk ''{print $0}'' %s',EVT);
[log, EVTSTR] = system(com);

%% Get station information
com = sprintf('awk ''{print $0}'' %s',STA);
[log, STASTR] = system(com);

%% Write file

if exist(RUNFILE,'file') == 2
    disp('File exists! Removing it now')
    com = ['rm -f',RUNFILE];
    [status,log] = system(com);
end

fid = fopen(RUNFILE,'w');
fprintf(fid,'%s\n',TABLE);               % path to mode table
fprintf(fid,'%.4f %.1f\n',LENGTH_HR,DT); % length (hours); 1/samplerate (s)
fprintf(fid,'35.\n');                    % reference frequency wref
fprintf(fid,'%s\n',CARD);                % card
fprintf(fid,'%s\n',COMP);                % T (Transverse; toroidal) Z (Vertical; spheroidal) R (Radial; spheroidal)
fprintf(fid,'.true.\n');                 % .true. = Displacement; .false. = Acceleration
fprintf(fid,'.false.\n');                % use crust5.1 model? May or may not work...
fprintf(fid,'sa_junk\n');                %junk if no crustal model
fprintf(fid,'frechet_file_junk\n');      %junk if no crustal model
fprintf(fid,'%s',EVTSTR);                % event file information
fprintf(fid,'%s',STASTR);                % station file information
fclose(fid);





    