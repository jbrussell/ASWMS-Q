% Invert beta gradient field to obtain the best fitting scalar beta field.
%
% This function solves equations 10-11 in Bao et al. (2016) GJI; doi:10.1093/gji/ggw151
% using least squares, with the additional constraint that the mean value
% of dlnbeta across the region is zero. This results in an average beta
% value of 1. An additional second derivative smoothing constraint is
% enforced
%
% % INPUT:
% xi = meshgrid (latitude)
% yi = mesghrid (longitude)
% dlnbetaLat_map = gradient in beta along x (lat) direction
% dlnbetaLon_map = gradient in beta along y (lon) direction
% smweight0 = smoothing weight (as fraction of G matrix norm)
%
% OUTPUT: 
% beta = desired scalar field values
%
% jbrussell 3/1/2022
function [beta_map]=inv_beta(xi,yi,dlnbetaLat_map,dlnbetaLon_map,smweight0)

[Nla,Nlo]=size(xi);
lnbeta=zeros(Nla,Nlo);
N = Nla*Nlo;

% Do latitudes first (x)
dlnbetaLat = zeros(N,1);
Glat = zeros(N,N);
for ila = 1:Nla
    for ilo = 1:Nlo
        ii=Nlo*(ila-1)+ilo;
        
        % Build data kernel
        if ila-1 == 0
            % endpoint, left derivative
            dla=vdist(xi(ila,ilo),yi(ila,ilo),xi(ila+1,ilo),yi(ila+1,ilo))/1e3;
            Glat(ii,ii) = -1 ./ dla;
            Glat(ii,ii+Nla) = 1 ./ dla;
        elseif ila+1 > Nla
            % endpoint, right derivative
            dla=vdist(xi(ila-1,ilo),yi(ila-1,ilo),xi(ila,ilo),yi(ila,ilo))/1e3;
            Glat(ii,ii-Nla) = -1 ./ dla;
            Glat(ii,ii) = 1 ./ dla;
        else
            dla1=vdist(xi(ila-1,ilo),yi(ila-1,ilo),xi(ila,ilo),yi(ila,ilo))/1e3;
            dla2=vdist(xi(ila,ilo),yi(ila,ilo),xi(ila+1,ilo),yi(ila+1,ilo))/1e3;
            Glat(ii,ii-Nla) = 1 ./ (dla1+dla2);
            Glat(ii,ii+Nla) = 1 ./ (dla1+dla2);
            Glat(ii,ii) = -2 ./ (dla1+dla2);
        end
        
        % Data vector
        dlnbetaLat(ii,1) = dlnbetaLat_map(ila,ilo);
    end
end

% Now longitudes (y)
dlnbetaLon = zeros(N,1);
Glon = zeros(N,N);
for ila = 1:Nla
    for ilo = 1:Nlo
        ii=Nlo*(ila-1)+ilo;
        
        % Build data kernel
        if ilo-1 == 0
            % endpoint, left derivative
            dlo=vdist(xi(ila,ilo),yi(ila,ilo),xi(ila,ilo+1),yi(ila,ilo+1))/1e3;
            Glon(ii,ii) = -1 ./ dlo;
            Glon(ii,ii+1) = 1 ./ dlo;
        elseif ilo+1 > Nlo
            % endpoint, right derivative
            dlo=vdist(xi(ila,ilo-1),yi(ila,ilo-1),xi(ila,ilo),yi(ila,ilo))/1e3;
            Glon(ii,ii-1) = -1 ./ dlo;
            Glon(ii,ii) = 1 ./ dlo;
        else
            dlo1=vdist(xi(ila,ilo-1),yi(ila,ilo-1),xi(ila,ilo),yi(ila,ilo))/1e3;
            dlo2=vdist(xi(ila,ilo),yi(ila,ilo),xi(ila,ilo+1),yi(ila,ilo+1))/1e3;
            Glon(ii,ii-1) = 1 ./ (dlo1+dlo2);
            Glon(ii,ii+1) = 1 ./ (dlo1+dlo2);
            Glon(ii,ii) = -2 ./ (dlo1+dlo2);
        end
        
        % Data vector
        dlnbetaLon(ii,1) = dlnbetaLon_map(ila,ilo);
    end
end


% Combine lat and lon
G = [Glat; Glon];
dlnbeta = [dlnbetaLat; dlnbetaLon];
% G = [Glat; Glon; Gconstraint];
% dlnbeta = [dlnbetaLat; dlnbetaLon; dconstraint];

% remove nan values;
inan = find(isnan(dlnbeta));
% igood = find(~isnan(dlnbeta));
G(inan,:) = [];
dlnbeta(inan,:) = [];
% % Remove unsampled grid points
% inode_nodata = find(sum(G,1)==0);
% inode_good = find(sum(G,1)~=0);
% G(:,inode_nodata) = [];

% Add extra constraint that the mean value equals 0
Nsolve = size(G,2);
Gconstraint = ones(1,Nsolve) ./ Nsolve;
dconstraint = 0;
G = [G; Gconstraint];
dlnbeta = [dlnbeta; dconstraint];

% Add smoothing kernel
xnode = xi(:,1)';
ynode = yi(1,:);
F = smooth_kernel_build(xnode, ynode, N);
f = zeros(size(F,1),1);
% Rescale the smooth kernel
NR=norm(F,1);
NA=norm(G,1);
smweight = smweight0*NA/NR;

% Construct full G amtrix
G = [G; smweight*F];
dlnbeta = [dlnbeta; smweight*f];

% Invert for receiver amplitude terms
% std_err = std(sqrt(dlnbetaLat_map(:).^2 + dlnbetaLon_map(:).^2)/5);
% W = diag(1./std_err).^2;
% F = W.^(0.5)*G;
% f = W.^(0.5)*dlnbeta;
% lnbeta = (F'*F)\F'*f;
% lnbeta = (G'*G)\G'*dlnbeta;
lnbeta = (G'*G)\G'*dlnbeta;

% Form matrix
beta_map = nan(size(xi));
for ila = 1:Nla
    for ilo = 1:Nlo
        ii=Nlo*(ila-1)+ilo;
        beta_map(ila,ilo) = exp(lnbeta(ii));
    end
end

if 0
    figure(999); clf
    ax = worldmap([min(xi(:)) max(xi(:))], [min(yi(:)) max(yi(:))]);
    surfacem(xi,yi,beta_map); colorbar;
    
    figure(1000); clf
    subplot(1,3,1);
    worldmap([min(xi(:)) max(xi(:))], [min(yi(:)) max(yi(:))]);
    surfacem(xi,yi,dlnbetaLat_map); colorbar;
    
    subplot(1,3,2);
    worldmap([min(xi(:)) max(xi(:))], [min(yi(:)) max(yi(:))]);
    surfacem(xi,yi,dlnbetaLon_map); colorbar;
    
    subplot(1,3,3);
    worldmap([min(xi(:)) max(xi(:))], [min(yi(:)) max(yi(:))]);
    surfacem(xi,yi,sqrt(dlnbetaLat_map.^2+dlnbetaLon_map.^2)); colorbar;
end

end

