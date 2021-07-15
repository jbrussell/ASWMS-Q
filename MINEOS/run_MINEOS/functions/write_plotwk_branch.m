%% Write driver for plot_wk
% NJA, 2014
% updated version in order to run for 200mHZ. 
% pylin.patty 2014/12,2015/01

function write_plotwk_branch(TYPE,CARDID,BR)

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

    RUNFILE = 'run_plotwk_branch.t';
    TYPEID = param.TTYPEID;

elseif strcmp(TYPE,'S') == 1
    
    RUNFILE = 'run_plotwk_branch.s';
    TYPEID = param.STYPEID;
    
else
    
    disp('Type does not exist!')

end
TABLEHDR = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.table_hdr'];
TABLEHDR_BR = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.table_hdr.branch'];
TABLE_MASK = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.',num2str(BR),'.mask'];

% Write out instructions to run file
fid = fopen(RUNFILE,'w');
fprintf(fid,'table %s\n',TABLEHDR);
fprintf(fid,'search \n');
fprintf(fid,'1 %4.2f %4.2f \n',[Wmin Wmax]);
fprintf(fid,'99 0 0 \n');
fprintf(fid,'branch %s\n',TABLEHDR_BR);
fprintf(fid,'search \n');
fprintf(fid,'1 %4.2f %4.2f \n',[Wmin Wmax]);
fprintf(fid,'99 0 0 \n');
fprintf(fid,'output %s b\n',TABLE_MASK);
fprintf(fid,'%d \n',BR);
fprintf(fid,'\n');
fprintf(fid,'quit \n');
fclose(fid);

% plot_wk <<!
% table $TABLEPATH$TABLE.table_hdr
% search
% 1 0. 201.
% 99 0 0
% branch $TABLEPATH$TABLE.table_hdr.branch
% search
% 1 0. 201.
% 99 0 0
% output $TABLEPATH$TABLE.$BR.mask b
% $BR
% 
% quit
% !
