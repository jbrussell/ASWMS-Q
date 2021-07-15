function [ FRECH ] = frech_aniso2( TYPE,CARDID,branch,periods,eig )
% Calculate Frechet Kernels for Anisotropy Inversion
% (Montagner & Nataf 1986)
%
% JBR - 2/14/2017
%
% JBR - 2/20/2017 :: Version 2 uses equations from Jim's thesis appendix
%

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

%% Load Card Model
card = read_model_card([param.CARDPATH,param.CARD]);
r = card.rad;
vsv = card.vsv;
vsh = card.vsh;
vpv = card.vpv;
vph = card.vph;
rho = card.rho;
eta = card.eta;

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
                
                
                if strcmp(TYPE,'T') == 1
                    wl = eig.ll(i).nn(j+1).wl;
                    wp = eig.ll(i).nn(j+1).wp;
                    ll = eig.ll(i).nn(j+1).ll;
                    %rho = eig.rho;
                    
                    % Jim's Thesis (A22)
                    FRECH(iper).L = rho.*vsv.^2.*ll*(ll+1).*(wp-wl./r).^2;
                    FRECH(iper).N = rho.*vsh.^2.*(ll+2).*(ll+1).*ll.*(ll-1).*(wl./r).^2;
                    FRECH(iper).rho = rho.*wl.^2;
                    FRECH(iper).per = per;
                    
                elseif strcmp(TYPE,'S') == 1
                    u = eig.ll(i).nn(j+1).u;
                    up = eig.ll(i).nn(j+1).up;
                    v = eig.ll(i).nn(j+1).v;
                    vp = eig.ll(i).nn(j+1).vp;
                    ll = eig.ll(i).nn(j+1).ll;
                    %rho = eig.rho;
                    
                    % Jim's Thesis (A21)
                    FRECH(iper).A = rho.*vph.^2.*(2*u-ll*(ll+1)*v).^2./r.^2;
                    FRECH(iper).C = rho.*vpv.^2.*up.^2;
                    FRECH(iper).F = rho.*eta.*(vph.^2-2*vsv.^2).*2.*up.*(2*u-ll*(ll+1)*v)./r;
                    FRECH(iper).L = rho.*vsv.^2.*ll*(ll+1).*(vp+(u-v)./r).^2;
                    %FRECH(iper).N = rho.*vsh.^2(ll+2
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

