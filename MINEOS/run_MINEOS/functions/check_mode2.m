% 10/7/15 -- Josh Russell
%
% Program to check for skipped and problematic eigenfrequencies from mineos.
% 
% INPUT: LOGF - log file from a mineos_nohang run.
%
% OUTPUT: l_start - the mode directly following the last successfully
% calculate mode.
%

function l_start = check_mode2(LOGF,LOOP,L_START_PREV)

parameter_FRECHET;
TYPE = param.TYPE;
MODE = ch_mode;
if strcmp(TYPE,'S') == 1
    TYPEID = param.STYPEID; 
elseif strcmp(TYPE,'T') == 1
    TYPEID = param.TTYPEID;
end

ASC = [CARDTABLE,param.CARDID,'.',TYPEID,'_',num2str(LOOP),'.asc'];

% check for the escape characters 'i quit' in ascii file
com = ['awk ''{if ($2=="quit") print $2}'' ',LOGF];
[log1, check] = system(com);
check = check(1:end-1);

% check for increasing period instead of decreasing


if strcmp(check,'quit') || 

    %extract last successfully calculated mode l and w;
    com = ['awk ''{ if ($1==',num2str(MODE),' && $2=="',TYPE,'") print $0 }'' ',ASC,' > ',CARDTABLE,'temp1'];
    [log2, temp1] = system(com);
    
    com = ['awk ''{print $3, $5}'' ',CARDTABLE,'temp1'];
    [log3, lw_all] = system(com);
    lw_all = str2num(lw_all);
    if isempty(lw_all) %Last mineos run produced no eigenfrequency calculations. skip
        disp(['No eigenfrequencies calculated for l= ',num2str(L_START_PREV)])
        l_start = L_START_PREV + 1;
        return
    end
    ll_all = lw_all(:,1);
    ww_all = lw_all(:,2);
    
    I_l = diff(ll_all)>1;
    check_lgap = ll_all(I_l);
    I_w = diff(ww_all)>1;
    check_wgap = ww_all(I_w);
    
    if ~isempty(check_lgap)
        if ~isempty(check_wgap)
            if ll_all(I_w) < ll_all(I_l)
                lgap = ll_all(I_w);
                l_start = lgap(1)+1; % start on first skipped w
            else
                lgap = ll_all(I_l);
                l_start = lgap(1)+1; % start on first skipped l
            end
        else
            l_start = check_lgap(1)+1; % start on first skipped l
        end
    elseif ~isempty(check_wgap)
        lgap = ll_all(I_w);
        l_start = lgap(1)+1; % start on first skipped w
    else
        l_start = ll_all(end)+1;
    end
else
    l_start = NaN;
end
    