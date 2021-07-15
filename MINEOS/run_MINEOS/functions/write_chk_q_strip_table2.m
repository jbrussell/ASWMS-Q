% 10/7/15 -- Josh Russell
% 
% Function to write the input files for mineos_q, mineos_strip, and
% mineos_table for each of the *.eig_fix files
%
% tables made with certain card files prefer this version for the q driver
% input... I do not know why.
%

function write_chk_q_strip_table2(LOOP)

parameter_FRECHET;
TYPE = param.TYPE;
CARDID = param.CARDID;
CARDPATH = param.CARDPATH;
RUNPATH = param.RUNPATH;
TABLEPATH = param.TABLEPATH;
MODEPATH = param.MODEPATH;
CARD = param.CARD;

% % Mineos table parameters
% maxN = 18000; % Estimate of max number of modes
% minF = 0;
% maxF = 200.1 %%150.05; %50.05;
% minL = 0;



CAR = [CARDPATH,CARD];

%% Set up the names of the files (diff. for spheroidal & toroidal)

if strcmp(TYPE,'S') == 1
    
    TYPEID = param.STYPEID;
    EIG0 = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'_0.eig_fix'];

    
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
    
elseif strcmp(TYPE,'T') == 1
    
    TYPEID = param.TTYPEID;
    EIG0 = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'_0.eig_fix'];

    
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
    
else
    disp('Type not recognized! Must be S or T');
end


%% Now write the files


% Write q driver
fid=fopen(RUNFILE_q,'w');
fprintf(fid,'%s\n',QMOD);
fprintf(fid,'%s\n',QOUT);
fprintf(fid,'%s\n',EIG0);
% fprintf(fid,'%s\n','y');
for i = 1:LOOP
    EIG = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'_',num2str(i),'.eig_fix'];
    fprintf(fid,'%s\n',EIG);
end
fprintf(fid,' \n');
fprintf(fid,'\n');
fclose(fid);

% Write strip driver
fid=fopen(RUNFILE_strip,'w');
fprintf(fid,'%s\n',STRIP);
fprintf(fid,'%s\n',EIG0);
for i = 1:LOOP
    EIG = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'_',num2str(i),'.eig_fix'];
    fprintf(fid,'%s\n',EIG);
end
fprintf(fid,'\n');
fclose(fid);

% Write table driver
fid=fopen(RUNFILE_table,'w');
fprintf(fid,'%s\n',TABLE);
fprintf(fid,'%s\n',num2str(maxN));
fprintf(fid,'%f %f\n',[minF maxF]);
fprintf(fid,'%i %i\n',[minL maxL]);
fprintf(fid,'%s\n',QOUT);
fprintf(fid,'%s\n',STRIP);
fprintf(fid,'\n');
fclose(fid);