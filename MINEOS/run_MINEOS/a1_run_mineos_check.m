% 10/7/15 -- JBR
%
% Run Minos Bran (modified 10/7/15)
% This program is similar to 'run_mineos' but checks for missed
% eigenfrequencies and loops until mode table finishes calculating all
% frequencies from minF to maxF.
%
% NJA, 2014
% pylin.patty 2014.. 
% JBR 2015
%
clear

%% get pamameters information 
parameter_FRECHET;
CARD = param.CARD;

%% Remove existing table files -- JOSH 8/24/15 to allow using mutliple branches
com = ['rm ',CARDTABLE,'*',lower(param.TYPE),num2str(floor(minF)),'to',num2str(floor(maxF)),'*'];
[status,log] = system(com);
com = ['rm ',CARDTABLE,'log*'];
[status,log] = system(com);

%% Set path to fortran executables
% setpath_mineos;


%% Change environment variables to deal with gfortran
setenv('GFORTRAN_STDIN_UNIT', '5') 
setenv('GFORTRAN_STDOUT_UNIT', '6') 
setenv('GFORTRAN_STDERR_UNIT', '0')


%% Run Spheroidal Branches First
if SONLY 
    if exist([param.MODEPATH,param.SMODEIN],'file') == 0
        write_mode_in(0,nan);
    end
    
    num_loop = 0;
    ll = [];
    
    disp('Beginning to Calculate Spheroidal Modes')
    disp(sprintf('\n'))
    TYPE = 'S';
    
    % Write out run files -- be sure that paths are correct!
    write_mineos_drivers(TYPE,CARD);
    
    %mineos_nohang for s1-s5
    disp('Running mineos_nohang - This will take several minutes');
    disp('S1');
    
    tic
    LOG = [CARDTABLE,'log',num2str(num_loop)];
    com = ['cat run_nohang.s | mineos_nohang > ',LOG];
    [status,log] = system(com);
% % %     fprintf(fid,log);
    if status ~= 0     
        error( 'something is wrong at mineos_nohang')
    end
    toc
    
            TYPEID = param.STYPEID;
        com = ['cat ',CARDTABLE,param.CARDID,'.',TYPEID,'.asc > ',CARDTABLE,param.CARDID,'.',TYPEID,'_',num2str(num_loop),'.asc'];
        [status,log] = system(com);
        com = ['cat ',CARDTABLE,param.CARDID,'.',TYPEID,'.eig > ',CARDTABLE,param.CARDID,'.',TYPEID,'_',num2str(num_loop),'.eig'];
        [status,log] = system(com);
        l_start = check_mode(LOG,num_loop,0); % Check that all eigenfrequencies were calculated
        mode_chk = l_start;
    while ~isnan(mode_chk)
        num_loop = num_loop + 1;
        ll = [ll; l_start];
        system('rm run_nohang.s');
        write_mode_in(l_start,num_loop); % Build new mode.in file starting from last successful w,l
        write_chk_mineos_nohang(TYPE,CARD,num_loop);
        
        disp(['--- Rerunning mineos_nohang: LOOP ',num2str(num_loop),' ---']);
        disp(['Starting at l = ',num2str(l_start)]);
        tic
        LOG = [CARDTABLE,'log',num2str(num_loop)];
        com = ['cat run_nohang.s | mineos_nohang > ',LOG];
        [status,log] = system(com);
        
        if status ~= 0     
            error( 'something is wrong at mineos_nohang')
        end
        toc
        
        l_start_prev = l_start;
        l_start = check_mode(LOG,num_loop,l_start_prev); % Check that all eigenfrequencies were calcualted 
        
        mode_chk = l_start;
%        pause;
    end
    
    % eig_recover
    if ~isempty(ll)
        disp(['Running eig_recover for ',num2str(num_loop),' files'])
        for i = 1:num_loop
            disp(['file ',num2str(i),' ...']);
            write_eig_recov(i-1,ll(i)-1);
            com = ['cat run_eigrecov.s | eig_recover'];
            [status log] = system(com);
        end
        com = ['cat ',CARDTABLE,param.CARDID,'.',TYPEID,'_',num2str(i),'.eig > ',CARDTABLE,param.CARDID,'.',TYPEID,'_',num2str(i),'.eig_fix'];
        [status,log] = system(com);     
    
        
        % WRITE DRIVERS
        write_chk_q_strip_table(num_loop);
    end
    
    % mineos_qcorrectphv
    disp('Running mineos_qcorrectphv');
    com = ['cat run_q.s | mineos_qcorrectphv'];
    [status,log] = system(com);
    if status ~= 0     
        error( 'something is wrong at mineos_qcorrectphv')
    end
    
    % mineos_strip
    disp('Running mineos_strip');
    com = ['cat run_strip.s | mineos_strip'];
    [status,log] = system(com);
    if status ~= 0     
        error( 'something is wrong at mineos_strip')
    end
    % mineos_table
    disp('Running mineos_table');
    com = ['cat run_table.s | mineos_table > log_table'];
    [status,log] = system(com);
    if status ~= 0     
        error( 'something is wrong at mineos_table')
    end
    
    delete('run_nohang.s','run_q.s', 'run_strip.s', 'run_table.s', 'run_eigrecov.s')
end

%% Run toroidal branches next
if TONLY
    if exist([param.MODEPATH,param.TMODEIN],'file') == 0
        write_mode_in(0,nan);
    end
    num_loop = 0;
    ll = [];
    
    disp('Beginning to Calculate Toroidal Modes')
    disp(sprintf('\n'));
    TYPE = 'T';
    % Write out run files -- be sure that paths are correct!
    write_mineos_drivers(TYPE,CARD);
    

    % mineos_nohang for s1-s5
    disp('Running mineos_nohang');
    tic
    LOG = [CARDTABLE,'log',num2str(num_loop)];
    com = ['cat run_nohang.t | mineos_nohang > ',LOG];
    [status,log] = system(com);
% % %     fprintf(fid,log);
    if status ~= 0     
        error( 'something is wrong at mineos_nohang')
    end
        toc
        
        TYPEID = param.TTYPEID;
        com = ['cat ',CARDTABLE,param.CARDID,'.',TYPEID,'.asc > ',CARDTABLE,param.CARDID,'.',TYPEID,'_',num2str(num_loop),'.asc'];
        [status,log] = system(com);
        com = ['cat ',CARDTABLE,param.CARDID,'.',TYPEID,'.eig > ',CARDTABLE,param.CARDID,'.',TYPEID,'_',num2str(num_loop),'.eig'];
        [status,log] = system(com);
        l_start = check_mode(LOG,num_loop,0); % Check that all eigenfrequencies were calculated
        mode_chk = l_start;
    while ~isnan(mode_chk)

        num_loop = num_loop + 1;
        ll = [ll; l_start];
        system('rm run_nohang.t');
        write_mode_in(l_start,num_loop); % Build new mode.in file starting from last successful w,l
        write_chk_mineos_nohang(TYPE,CARD,num_loop);
        
        disp(['--- Rerunning mineos_nohang: LOOP ',num2str(num_loop),' ---']);
        disp(['Starting at l = ',num2str(l_start)]);
%         tic
        LOG = [CARDTABLE,'log',num2str(num_loop)];
        com = ['cat run_nohang.t | mineos_nohang > ',LOG];
        [status,log] = system(com);
        
        if status ~= 0     
            error( 'something is wrong at mineos_nohang')
        end
%         toc
        
        l_start_prev = l_start;
        l_start = check_mode(LOG,num_loop,l_start_prev); % Check that all eigenfrequencies were calcualted 
        
        mode_chk = l_start;
%        pause;
    end
    
    % eig_recover
    if ~isempty(ll)
        disp(['Running eig_recover for ',num2str(num_loop),' files'])
        for i = 1:num_loop
            disp(['file ',num2str(i),' ...']);
            write_eig_recov(i-1,ll(i)-1);
            com = ['cat run_eigrecov.t | eig_recover'];
            [status log] = system(com);
        end
        com = ['cat ',CARDTABLE,param.CARDID,'.',TYPEID,'_',num2str(i),'.eig > ',CARDTABLE,param.CARDID,'.',TYPEID,'_',num2str(i),'.eig_fix'];
        [status,log] = system(com);     
    
        
        % WRITE DRIVERS
        write_chk_q_strip_table(num_loop);
    end
    
    % mineos_qcorrectphv
    disp('Running mineos_qcorrectphv');
    com = ['cat run_q.t | mineos_qcorrectphv > qlog'];
    [status,log] = system(com);
    if status ~= 0     
        write_chk_q_strip_table2(num_loop); % Some cards with *.eig_fix files prefer this version...
        disp('Something broke... trying mineos_qcorrectphv again');
        com = ['cat run_q.t | mineos_qcorrectphv > qlog'];
        [status,log] = system(com);
        if status ~= 0
            error( 'something is wrong at mineos_qcorrectphv')
        end
    end
    
    % mineos_strip
    disp('Running mineos_strip');
    com = ['cat run_strip.t | mineos_strip'];
    [status,log] = system(com);
    if status ~= 0 
        error( 'something is wrong at mineos_strip')
    end
    % mineos_table
    disp('Running mineos_table');
    com = ['cat run_table.t | mineos_table > mineos_table.LOG'];
    [status,log] = system(com);
    if status ~= 0     
        error( 'something is wrong at mineos_table')
    end
    
    delete('run_nohang.t','run_q.t', 'run_strip.t', 'run_table.t', 'run_eigrecov.t')

end

% Delete unnecessary files
delete('*.LOG','qlog','log_table');
delete([CARDTABLE,'log*']);
delete([CARDTABLE,'*.asc']);
if exist([CARDTABLE,param.CARDID,'.',TYPEID,'_0.eig_fix'],'file') ~= 0
    delete([CARDTABLE,'*',TYPEID,'*.eig']); % Delete regular .eig files
elseif num_loop == 0
    delete([CARDTABLE,param.CARDID,'.',TYPEID,'_',num2str(num_loop),'.eig']);
end

%% Change the environment variables back to the way they were
setenv('GFORTRAN_STDIN_UNIT', '-1') 
setenv('GFORTRAN_STDOUT_UNIT', '-1') 
setenv('GFORTRAN_STDERR_UNIT', '-1')

