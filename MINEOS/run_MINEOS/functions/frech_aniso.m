function [ FRECH ] = frech_aniso( TYPE,CARDID,branch,periods,eig )
% Calculate Frechet Kernels for Anisotropy Inversion
% (Montagner & Nataf 1986)
%
% JBR - 2/14/2017

% Get useful info from parameter file
parameter_FRECHET;
CARDPATH = param.CARDPATH;
FRECHETPATH = param.frechetpath;
TABLEPATH = param.TABLEPATH;
%periods = param.periods;

if strcmp(TYPE,'T') == 1
    disp('Toroidal!');
    TYPEID = param.TTYPEID;    
elseif strcmp(TYPE,'S') == 1
    disp('Spheroidal!');
    TYPEID = param.STYPEID;    
else
    disp('No TYPE recognized!');
    
end

%% Load dispersion curve
QIN = [TABLEPATH,CARDID,'/tables/',CARDID,'.',TYPEID,'.q'];
dat = {};
for i = branch+1
    com = ['awk ''{ if ($1 ==',num2str(i-1),' && $10 != "") print $0}'' ',QIN];
    [log3, dat{i}] = system(com);
    dat{i} = str2num(dat{i});
    nn =  dat{i}(:,1);
    ll =  dat{i}(:,2);
    w =   dat{i}(:,3)/(2*pi)*1000; %convert rad/s ---> mhz
    qq =  dat{i}(:,4);
    phi = dat{i}(:,5);
    cv =  dat{i}(:,6);
    gv =  dat{i}(:,7);
    cvq = dat{i}(:,8);
    Tq =  dat{i}(:,9);
    T =   dat{i}(:,10);
end

%% Calculate Kernels
num_ll = length(eig.ll);
iper = 0;
check_pers = zeros(1,length(periods));
for i = 1:1:num_ll
    for j = branch
        if (size(eig.ll(i).nn,2) - 1) >= j
        if ~isempty(eig.ll(i).nn(j+1).ll) && ~isempty(eig.ll(i).nn)
            per = eig.ll(i).nn(j+1).per;

            Iper = find(abs(per-periods)<per*.005);
            if ~isempty(Iper)
            if abs(per-periods(Iper))<per*.005 && check_pers(Iper)==0
                display([num2str(periods(Iper)),' s --> ',num2str(per)])
                check_pers(Iper) = 1;
                iper = iper + 1;
                
                [~, Icvq] = min(abs(periods(Iper)-T));
                c = cvq(Icvq);
                
                k = 2*pi/per/c; % wavenumber
                
                if strcmp(TYPE,'T') == 1
                    wl = eig.ll(i).nn(j+1).wl;
                    wp = eig.ll(i).nn(j+1).wp;
                    rho = eig.rho;
                    
                    % (Montagner & Nataf eq. 2)
                    FRECH(iper).L = (wp.^2)/k^2;
                    FRECH(iper).N = wl.^2;
                    FRECH(iper).rho = rho.*wl.^2;
                    FRECH(iper).per = per;
                    
                elseif strcmp(TYPE,'S') == 1
                    u = eig.ll(i).nn(j+1).u;
                    up = eig.ll(i).nn(j+1).up;
                    v = eig.ll(i).nn(j+1).v;
                    vp = eig.ll(i).nn(j+1).vp;
                    rho = eig.rho;
                    
                    % (Montagner & Nataf eq. 4)
                    FRECH(iper).A = v.^2;
                    FRECH(iper).C = up.^2/k^2;
                    FRECH(iper).F = 2*up.*v/k;
                    FRECH(iper).L = (vp/k-u).^2;
                    FRECH(iper).rho = rho.*(u.^2+v.^2);
                    FRECH(iper).per = per;
                end
                
            end
            end
        end
        end
    end
end

end

