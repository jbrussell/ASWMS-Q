function [mat] = load_azianifile(file,xnode,ynode)
%Load azimuthal anisotropy file and interpolate to desired grid
% Anisotropy amplitude output from TPW is zero-to-peak in km/s
%
isplot = 0;

[xi_want, yi_want] = ndgrid(xnode,ynode);

data = load(file);

Ac = data(:,1); % km/s
Ac_std = data(:,2); % km/s
As = data(:,3); % km/s
As_std = data(:,4); % km/s
lon = data(:,5);
lat = data(:,6);

if length(lat) > 1
    [Ac_mat] = griddata(lat,lon,Ac,xi_want,yi_want);
    [Ac_std_mat] = griddata(lat,lon,Ac_std,xi_want,yi_want);
    [As_mat] = griddata(lat,lon,As,xi_want,yi_want);
    [As_std_mat] = griddata(lat,lon,As_std,xi_want,yi_want);
else
    xi_want = lat;
    yi_want = lon;
    Ac_mat = Ac;
    Ac_std_mat = Ac_std;
    As_mat = As;
    As_std_mat = As_std;
end

% Calculate anisotropy parameters
A2 = sqrt(Ac_mat.^2 + As_mat.^2);
phi2 = 0.5*atan2d(As_mat,Ac_mat);

% Propagate error
A2_std = sqrt( ((Ac_mat.*Ac_std_mat).^2 + (As_mat.*As_std_mat).^2) ./ (Ac_mat.^2 + As_mat.^2) );
phi2_std = 180/pi * 0.5*sqrt( (As_mat.*Ac_std_mat).^2 + (Ac_mat.*As_std_mat).^2 ) ./ (Ac_mat.^2 + As_mat.^2);

mat.A2_kms = A2; % km/s
mat.A2_std_kms = A2_std; % km/s
mat.phi2 = phi2;
mat.phi2_std = phi2_std;
mat.Ac_kms = Ac_mat; % km/s
mat.Ac_std_kms = Ac_std_mat; % km/s
mat.As_kms = As_mat; % km/s
mat.As_std_kms = As_std_mat; % km/s
mat.lat = xi_want;
mat.lon = yi_want;

if isplot
    figure(200); clf;
    surface(xi_want-(xnode(2)-xnode(1))/2,yi_want-(ynode(2)-ynode(1))/2,zeros(size(xi_want)),A2); hold on;
end

end

