%% Write Driver Files for mineos_nohang,mineos_q,mineos_strip,mineos_table
% NJA, 2014

% There will be 7 of these (6 for radial and 1 for toroidal).
% TYPE is either S (spheroidal) or T (toroidal)
%
% 8/2014 Modified to use only one run file for spheroidal
%
% 10/27/2014 Modified to write out run files for ambient noise period bands
% Toroidal modes can be done with a single run file while spheroidal modes
% have to be split into two.  
%
% Allowed types are 'T' and 'S' -- S will create 2 run files
% NJA, 2014
% pylin.patty 2014
% updated version in order to run for 200mHZ. 
% pylin.patty 2015/01


function write_mineos_drivers(TYPE,CARD)

parameter_FRECHET;
CARDID = param.CARDID;
CARDPATH = param.CARDPATH;
RUNPATH = param.RUNPATH;
TABLEPATH = param.TABLEPATH;
MODEPATH = param.MODEPATH;

% % Mineos table parameters
% maxN = 18000; % Estimate of max number of modes
% minF = 0;
% maxF = 200.1 %%150.05; %50.05;
% minL = 0;



CAR = [CARDPATH,CARD];

%% Set up the names of the files (diff. for spheroidal & toroidal)

if strcmp(TYPE,'S') == 1
    
    %TYPEID = 's0to200';%'s0to50';
    %MODENUM = 's.mode0_200';%'s.mode0_50';
    TYPEID = param.STYPEID;
    MODENUM = param.SMODEIN;
    
    % Nohang File
    RUNFILE_nohang = 'run_nohang.s';

    
    % Name output files
    ASC = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.asc'];
    EIG = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.eig'];
    MOD = [MODEPATH,MODENUM];
    
    %% mineos_q driver
    RUNFILE_q = 'run_q.s';

    QMOD = [CARDPATH,CARDID,'.qmod'];
    QOUT = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.q'];

    %% Mineos Strip Driver
    RUNFILE_strip ='run_strip.s';
    STRIP = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.strip'];

    %% Mineos Table Driver
    RUNFILE_table = 'run_table.s';
    TABLE = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.table'];

    %% See mineos_table.f for more information about following parameters
%    maxL = 6000; %680;
    
%    %Write out short period S run file
%     RUNFILE_nohang2 = 'run_nohang2.s';
%     TYPEID2 = 's50to150';
%     MODENUM2 = 's.mode50_150';
%     ASC2 = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID2,'.asc'];
%     EIG2 = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID2,'.eig'];
%     MOD2 = [MODEPATH,MODENUM2];
% 
%     fid=fopen(RUNFILE_nohang2,'w');
%     fprintf(fid,'%s\n',CAR);
%     fprintf(fid,'%s\n',ASC2);
%     fprintf(fid,'%s\n',EIG2);
%     fprintf(fid,'%s\n',MOD2);
%     fprintf(fid,'\n');
%     fclose(fid);

    
elseif strcmp(TYPE,'T') == 1
    %TYPEID = 't0to150';
    %MODENUM = 't.mode0_150';
    
    TYPEID = param.TTYPEID;
    MODENUM = param.TMODEIN;
    % Nohang File
    RUNFILE_nohang = 'run_nohang.t';;

    
    % nohang driver
    ASC = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.asc'];
    EIG = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.eig'];
    MOD = [MODEPATH,MODENUM];
    
    
        % Name output files
    ASC = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.asc'];
    EIG = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.eig'];
    MOD = [MODEPATH,MODENUM];
    
    %% mineos_q driver
    RUNFILE_q = 'run_q.t';

    QMOD = [CARDPATH,CARDID,'.qmod'];
    QOUT = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.q'];

    %% Mineos Strip Driver
    RUNFILE_strip ='run_strip.t';
    STRIP = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.strip'];

    %% Mineos Table Driver
    RUNFILE_table = 'run_table.t';
    TABLE = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.table'];

    %% See mineos_table.f for more information about following parameters
%    maxL = 6000;%600; 


    
else
    disp('Type not recognized! Must be S or T');
end

%% Now write the files

% Write nohang driver
fid=fopen(RUNFILE_nohang,'w');
fprintf(fid,'%s\n',CAR);
fprintf(fid,'%s\n',ASC);
fprintf(fid,'%s\n',EIG);
fprintf(fid,'%s\n',MOD);
fprintf(fid,'\n');
fclose(fid);


% Wrtie q driver
fid=fopen(RUNFILE_q,'w');
fprintf(fid,'%s\n',QMOD);
fprintf(fid,'%s\n',QOUT);
fprintf(fid,'%s\n',EIG);
fprintf(fid,'%s\n','y');
%fprintf(fid,'%s\n',EIG2);
fprintf(fid,'\n');
fclose(fid);

% Write strip driver
fid=fopen(RUNFILE_strip,'w');
fprintf(fid,'%s\n',STRIP);
fprintf(fid,'%s\n',EIG);
%fprintf(fid,'%s\n',EIG2);
fprintf(fid,'\n');
fclose(fid);

% Write table driver
fid=fopen(RUNFILE_table,'w');
fprintf(fid,'%s\n',TABLE);
fprintf(fid,'%s\n',num2str(maxN));
fprintf(fid,'%i %f\n',[minF maxF]);
fprintf(fid,'%i %i\n',[minL maxL]);
fprintf(fid,'%s\n',QOUT);
fprintf(fid,'%s\n',STRIP);
fprintf(fid,'\n');
fclose(fid);

% Now move the files to the run directory
com = ['mv run_* ',RUNPATH,'.'];
[status,log] = system(com);
