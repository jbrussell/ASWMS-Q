% load_frechcv_asc.
% loads frechet kernels from ascii files
%
% JBR 10/11/16
%

function [FRECH] = load_frechcv_asc_per(TYPE,CARD,BRANCH,periods)

% Get useful info from parameter file
parameter_FRECHET;
CARDPATH = param.CARDPATH;
FRECHETPATH = param.frechetpath;
TABLEPATH = param.TABLEPATH;
% periods = param.periods;

if strcmp(TYPE,'T') == 1
    disp('Toroidal!');
    
    TYPEID = param.TTYPEID;
    
elseif strcmp(TYPE,'S') == 1
    disp('Spheroidal!');
    
    TYPEID = param.STYPEID;
    
else
    disp('No TYPE recognized!');
    
end

BRID = [num2str(BRANCH)];

for ip = 1:length(periods)
    
    FRECHASC = [FRECHETPATH,CARD,'.',TYPEID,'.',BRID,'.',num2str(periods(ip))];


    
    % Load in frechet files for each period
    
    % spheroidal, no aniso: 1=Vs,2=Vp,3=rho
    % spheroidal, aniso: 1=Vsv,2=Vpv,3=Vsh,4=Vph,5=eta,6=rho
    % toroidal, no aniso: 1=Vs,2=rho
    % toroidal, aniso: 1=Vsv,2=Vsh,3=rho
    % disp(FRECHASC)
    fid = fopen(FRECHASC,'r');
    C = textscan(fid,'%f%f%f%f%f%f%f');
    
    if strcmp(TYPE,'S') == 1
        FRECH(ip).per = periods(ip);
        FRECH(ip).rad = C{1};
        FRECH(ip).vsv = C{2};
        FRECH(ip).vpv = C{3};
        FRECH(ip).vsh = C{4};
        FRECH(ip).vph = C{5};
        FRECH(ip).eta = C{6};
        FRECH(ip).rho = C{7};
    elseif strcmp(TYPE,'T') == 1
        FRECH(ip).per = periods(ip);
        FRECH(ip).rad = C{1};
        FRECH(ip).vsv = C{2};
        FRECH(ip).vsh = C{3};
        FRECH(ip).rho = C{4};
    end
end


end