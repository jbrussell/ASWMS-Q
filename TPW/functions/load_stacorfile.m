function [mat] = load_stacorfile(file,stamat)
%Load alpha file and interpolate to desired grid
%
data = readtable(file);

sta = data.Var4;
stacor = data.Var2;
stacor_std = data.Var3;
staid = data.Var1;

stas = load(stamat);
[~,I] = intersect(stas.stnms,sta);
stlas = stas.stlas(I);
stlos = stas.stlos(I);

mat.sta = sta;
mat.stacor = stacor;
mat.stacor_std = stacor_std;
mat.staid = staid;
mat.stlas = stlas;
mat.stlos = stlos;

end

