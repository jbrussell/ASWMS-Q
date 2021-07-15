%% Write driver for plot_wk
% NJA, 2014
% updated version in order to run for 200mHZ. 
% pylin.patty 2014/12,2015/01

function write_plotwk(TYPE,CARDID)

parameter_FRECHET;
TABLEPATH = param.TABLEPATH;
% Parameters for searching for modes
%Wmin = 0;
%Wmax = 200;%150;%51;
Wmin = minF;
Wmax = maxF;

% TYPE = 'T';
% CARD = param.CARDID;

if strcmp(TYPE,'T') == 1

    RUNFILE = 'run_plotwk.t';
    TYPEID = param.TTYPEID;
    TABLEHDR = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.table_hdr'];
    %TABLEHDR = [TABLEPATH,CARDID,'/tables/',CARDID,'.t0to150.tab_hdr'];

elseif strcmp(TYPE,'S') == 1
    
    RUNFILE = 'run_plotwk.s';
    TYPEID = param.STYPEID;
    TABLEHDR = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.table_hdr'];
    %TABLEHDR = [TABLEPATH,CARDID,'/tables/',CARDID,'.s0to50.tab_hdr'];
    
else
    
    disp('Type does not exist!')

end

% Write out instructions to run file
fid = fopen(RUNFILE,'w');
fprintf(fid,'table %s\n',TABLEHDR);
fprintf(fid,'search \n');
fprintf(fid,'1 %4.2f %4.2f \n',[Wmin Wmax]);
fprintf(fid,'99 0 0 \n');
fprintf(fid,'branch \n');
fprintf(fid,'\n');
fprintf(fid,'quit \n');
fclose(fid);
