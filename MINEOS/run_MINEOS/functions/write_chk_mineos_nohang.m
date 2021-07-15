% 10/7/15 -- Josh Russell
%
% Write Driver File for mineos_nohang after checking for missing eigenfrequencies
%


function write_chk_mineos_nohang(TYPE,CARD,LOOP)

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
    ASC = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'_',num2str(LOOP),'.asc'];
    EIG = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'_',num2str(LOOP),'.eig'];
    MOD = [MODEPATH,MODENUM,'_',num2str(LOOP)];
    

    
elseif strcmp(TYPE,'T') == 1
    %TYPEID = 't0to150';
    %MODENUM = 't.mode0_150';
    
    TYPEID = param.TTYPEID;
    MODENUM = param.TMODEIN;
    % Nohang File
    RUNFILE_nohang = 'run_nohang.t';
    
    
        % Name output files
    ASC = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'_',num2str(LOOP),'.asc'];
    EIG = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'_',num2str(LOOP),'.eig'];
    MOD = [MODEPATH,MODENUM,'_',num2str(LOOP)];

    
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