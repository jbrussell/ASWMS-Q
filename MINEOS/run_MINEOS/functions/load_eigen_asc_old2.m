% Written on 9/14/15: Josh Russell
%
% Function converts raw (unformatted) eigenfunction file output from Mineos
% into an ascii file using forgran code 'eigenST_asc'.
%
% MUST RUN mk_kernels FIRST !!!!!
% 
% OUTPUT DESCRIPTIONS: 
%
%        eigen -- full data matrix including all angular orders (ll) and
%        modes (nn)
%
% SPHEROIDAL:
%         angular oder ->       ll = eigen(:,1)
%         mode number  ->       nn = eigen(:,2)
%         angular freq ->       w = eigen(:,3)
%         radius (m)    ->      r = eigen(:,4)
%         bulk modulus kappa -> kappa = eigen(:,5)
%         shear modulus mu ->   mu = eigen(:,6)
%         density rho   ->      dnn = eigen(:,7)
%         vert eigenf  ->       u = eigen(:,8);
%         d/dr(u)         ->    up = eigen(:,9);
%         horiz eigenf ->       v = eigen(:,10);
%         d/dr(v)         ->    vp = eigen(:,11);
%         ??? eigenf   ->       phi = eigen(:,12);
%         d/dr(phi)        ->   phip = eigen(:,13);
%         vertical stress ->    R
%         horiz stress ->       S 
%         radius to plot   ->   rad = 6371 - eigen(:,4)
%         period   ->           per = 2*pi/w
%
% TOROIDAL:
%         angular oder ->       ll = eigen(:,1)
%         mode number  ->       nn = eigen(:,2)
%         angular freq ->       w = eigen(:,3)
%         radius (m)    ->      r = eigen(:,4)
%         bulk modulus kappa -> kappa = eigen(:,5)
%         shear modulus mu ->   mu = eigen(:,6)
%         density rho   ->      dnn = eigen(:,7)
%         horiz eigenf ->       wl = eigen(:,8);
%         d/dr(wl)         ->   wp = eigen(:,9);
%         horiz stress ->       T 
%         radius to plot   ->   rad = 6371 - eigen(:,4)
%         period   ->           per = 2*pi/w
%
%        eig -- data structure in the form eig.ll(i).nn(j).dat where i is
%        the angular order and j is 1 + mode #. For example
%        eig.ll(10).nn(1).dat is the 10th angular order for the 0th
%        (fundamental) mode.
%
% JBR 10/6/16 -- Modified to handle *.eig_fix files
%
% JBR 1/20/17 -- Modified to include both Spheroidal and Toroidal
%

function [eigen,eig,saveopt] = load_eigen_asc(TYPE)

% Get useful info from parameter file
parameter_FRECHET;
CARDPATH = param.CARDPATH;
TABLEPATH = param.TABLEPATH;
CARDID = param.CARDID;
if ( TYPE == 'T') 
    TYPEID = param.TTYPEID;
elseif ( TYPE == 'S') 
    TYPEID = param.STYPEID;
end
EIGPATH = param.eigpath;

setpath_mineos;


% Normalization for eigenfunctions
% scale  = 1.0/(rn*sqrt(rn*pi*bigg)*rhobar)
rn = 6371000;
bigg = 6.6723e-11; % m/kg/s2
rhobar = 5515; % kg/m3
scale = 1/(rn*sqrt(rn*pi*bigg)*rhobar);
scale = 1;

%% Change environment variables to deal with gfortran
setenv('GFORTRAN_STDIN_UNIT', '5') 
setenv('GFORTRAN_STDOUT_UNIT', '6') 
setenv('GFORTRAN_STDERR_UNIT', '0')

%% Setup runfile
NDISC = 0;
if strcmp(TYPE,'T') == 1
    disp('Toroidal!');
    BRID = 2;
elseif strcmp(TYPE,'S') == 1
    disp('Spheroidal!');
    BRID = 3;
end

BRANCH = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.table_hdr.branch'];
EIG_ASC = [EIGPATH,CARDID,'.',TYPEID,'.',num2str(N_modes),'.asc'];
if exist(EIG_ASC,'file') == 2
%disp('File exists! Removing it now')
com = ['rm -f ',EIG_ASC];
[status,log] = system(com);
end

eig_fils = dir([param.TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'*.eig_fix']);
nfils = size(eig_fils,1);

if nfils == 0 % %no .eig_fix files
    EIG_I = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.eig'];

    RUNFILE = [EIGPATH,'eigen.in'];
    % Write runfile for eigen_asc
    fid = fopen(RUNFILE,'w');
    fprintf(fid,'%s\n',BRANCH);%input branch file
    fprintf(fid,'%s\n',EIG_ASC); %output ascii file
    fprintf(fid,'%s\n',EIG_I); %input raw .eig file
    fprintf(fid,'%d\n',NDISC);
    fclose(fid);
    
    com = sprintf('cat %s | eigenST_asc',RUNFILE);
    [status,log] = system(com);
    
elseif nfils ~= 0    
% LOOP OVER *.eig_fix FILES
for ifil = 1:nfils
    EIG_ASC_FIX = [EIGPATH,CARDID,'.',TYPEID,'.',num2str(N_modes),'_',num2str(ifil-1),'_fix.asc'];
    EIG_I = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'_',num2str(ifil-1),'.eig_fix'];
    display(['Working on: ',EIG_I])

    RUNFILE = [EIGPATH,'eigen.in'];
    if exist(RUNFILE,'file') == 2
    %disp('File exists! Removing it now')
    com = ['rm -f ',RUNFILE];
    [status,log] = system(com);
    end
    % Write runfile for eigen_asc
    fid = fopen(RUNFILE,'w');
    fprintf(fid,'%s\n',BRANCH);%input branch file
    fprintf(fid,'%s\n',EIG_ASC_FIX); %output ascii file
    fprintf(fid,'%s\n',EIG_I); %input raw .eig file
    fprintf(fid,'%d\n',NDISC);
    fclose(fid);

    %% Run fortran program

    if exist(EIG_ASC_FIX,'file') == 2
    %disp('File exists! Removing it now')
    com = ['rm -f ',EIG_ASC_FIX];
    [status,log] = system(com);
    end

    com = sprintf('cat %s | eigenST_asc',RUNFILE);
%     com = sprintf('cat %s | eigenSTscaled_asc',RUNFILE);
    [status,log] = system(com);
    
    % COMBINE ASCI FILES
    system(['cat ',EIG_ASC_FIX,' >> ',EIG_ASC]);
    delete(EIG_ASC_FIX);
end
end

%% Load eigenfunctions into data structure

% SPHEROIDAL
if (TYPE == 'S')
    eigen = load(sprintf('%s',EIG_ASC));
    ll = eigen(:,1);
    nn = eigen(:,2);
    w = eigen(:,3);
    r = eigen(:,4);
    kappa = eigen(:,5);
    mu = eigen(:,6);
    rho = eigen(:,7);
    u = eigen(:,8)*scale;
    up = eigen(:,9)*scale;
    v = eigen(:,10)*scale;
    vp = eigen(:,11)*scale;
    phi = eigen(:,12);
    phip = eigen(:,13);
    
    % CALCULATE vert traction R (Dahlen & Tromp eq. 8.47) k = sqrt(l*(l+1)
    R = (kappa+4/3*mu).*up + (kappa-2/3*mu).*(2*u-sqrt(ll.*(ll+1)).*v)./r;
    % Horizontal traction S (D&T eq. 8.48)
    S = mu.*(vp-v./r + sqrt(ll.*(ll+1)).*u./r);
    
    eigen(:,14) = R(:,:);
    eigen(:,15) = S(:,:);
    
    % Sort eigenfunctions into structure
    num_ll = ll(end);
    for i = 1:num_ll %loop over angular orders
        eigen_ll = eigen(eigen(:,1) == i,:);
        if ~isempty(eigen_ll)
            first_nn = eigen_ll(1,2);
            last_nn = eigen_ll(end,2);
            for j = first_nn:last_nn %loop over modes
                eigen_ll_nn = eigen_ll(eigen_ll(:,2) == j,:);
                %eig.nn(j+1).ll(i).dat = eigen_ll_nn(:,:);
                eig.nn(j+1).ll(i).ll = eigen_ll_nn(1,1);
                eig.nn(j+1).ll(i).nn = eigen_ll_nn(1,2);
                eig.nn(j+1).ll(i).w = eigen_ll_nn(1,3);
                eig.r = eigen_ll_nn(:,4);
                eig.kappa = eigen_ll_nn(:,5);
                eig.mu = eigen_ll_nn(:,6);
                eig.rho = eigen_ll_nn(:,7);
                eig.nn(j+1).ll(i).u = eigen_ll_nn(:,8)*scale;
                eig.nn(j+1).ll(i).up = eigen_ll_nn(:,9)*scale;
                eig.nn(j+1).ll(i).v = eigen_ll_nn(:,10)*scale;
                eig.nn(j+1).ll(i).vp = eigen_ll_nn(:,11)*scale;
                eig.nn(j+1).ll(i).phi = eigen_ll_nn(:,12);
                eig.nn(j+1).ll(i).phip = eigen_ll_nn(:,13);
                eig.nn(j+1).ll(i).R = eigen_ll_nn(:,14);
                eig.nn(j+1).ll(i).S = eigen_ll_nn(:,15);
                eig.rad = 6371 - eig.r/1000;
                eig.nn(j+1).ll(i).per = 2*pi./eig.nn(j+1).ll(i).w;
            end
        else
            eig.ll(i).nn = [];
        end
    end
    
% TOROIDAL
elseif ( TYPE == 'T') 
    eigen = load(sprintf('%s',EIG_ASC));
    ll = eigen(:,1);
    nn = eigen(:,2);
    w = eigen(:,3);
    r = eigen(:,4);
    kappa = eigen(:,5);
    mu = eigen(:,6);
    rho = eigen(:,7);
    wl = eigen(:,8)*scale;
    wp = eigen(:,9)*scale;
    
    % CALCULATE horiz traction T (Dahlen & Tromp eq. 8.49)
    T = mu.*(wp-wl./r);
    
    eigen(:,10) = T(:,:);
    
    % Sort eigenfunctions into structure
    num_ll = ll(end);
    for i = 1:num_ll %loop over angular orders
        eigen_ll = eigen(eigen(:,1) == i,:);
        if ~isempty(eigen_ll)
            first_nn = eigen_ll(1,2);
            last_nn = eigen_ll(end,2);
            for j = first_nn:last_nn %loop over modes
                eigen_ll_nn = eigen_ll(eigen_ll(:,2) == j,:);
                %eig.nn(j+1).ll(i).dat = eigen_ll_nn(:,:);
                eig.nn(j+1).ll(i).ll = eigen_ll_nn(1,1);
                eig.nn(j+1).ll(i).nn = eigen_ll_nn(1,2);
                eig.nn(j+1).ll(i).w = eigen_ll_nn(1,3);
                eig.r = eigen_ll_nn(:,4);
                eig.kappa = eigen_ll_nn(:,5);
                eig.mu = eigen_ll_nn(:,6);
                eig.rho = eigen_ll_nn(:,7);
                eig.nn(j+1).ll(i).wl = eigen_ll_nn(:,8)*scale;
                eig.nn(j+1).ll(i).wp = eigen_ll_nn(:,9)*scale;
                eig.nn(j+1).ll(i).T = eigen_ll_nn(:,10);
                eig.rad = 6371 - eig.r/1000;
                eig.nn(j+1).ll(i).per = 2*pi./eig.nn(j+1).ll(i).w;
            end
        else
            eig.ll(i).nn = [];
        end
    end
    
end

% Check if size of 'eig' is greater than 2 Gb
varinfo=whos('eig');
saveopt='';
if varinfo.bytes >= 2^31
  saveopt='-v7.3';
end

% Change the environment variables back to the way they were
setenv('GFORTRAN_STDIN_UNIT', '-1') 
setenv('GFORTRAN_STDOUT_UNIT', '-1') 
setenv('GFORTRAN_STDERR_UNIT', '-1')