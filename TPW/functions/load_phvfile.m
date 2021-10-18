function [mat] = load_phvfile(file,xnode,ynode)
%Load phase velocity file and interpolate to desired grid
%
isplot = 0;

[xi_want, yi_want] = ndgrid(xnode,ynode);

data = load(file);

phv = data(:,2);
phv_std = data(:,3);
resid_diag = data(:,4);
lon = data(:,5);
lat = data(:,6);

if length(lat) > 1
    [phv_mat] = griddata(lat,lon,phv,xi_want,yi_want);
    [phv_std_mat] = griddata(lat,lon,phv_std,xi_want,yi_want);
    [resid_diag_mat] = griddata(lat,lon,resid_diag,xi_want,yi_want);
else
    xi_want = lat;
    yi_want = lon;
    phv_mat = phv;
    phv_std_mat = phv_std;
    resid_diag_mat = resid_diag;
end

mat.phv = phv_mat;
mat.phv_std = phv_std_mat;
mat.resid_diag = resid_diag_mat;
mat.lat = xi_want;
mat.lon = yi_want;

if isplot
    figure(200); clf;
    surface(xi_want-(xnode(2)-xnode(1))/2,yi_want-(ynode(2)-ynode(1))/2,zeros(size(xi_want)),phv_mat); hold on;
    scatter(lat,lon,150,phv,'filled','MarkerEdgeColor','k')
end

end

