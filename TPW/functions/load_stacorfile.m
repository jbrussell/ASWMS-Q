function [mat] = load_stacorfile(file,stamat)
%Load alpha file and interpolate to desired grid
%
fid = fopen(file);
ista = 0;
while ~feof(fid)
    ista = ista + 1;
    tline = fgetl(fid);
    scan = textscan(tline,'%f %f %f %s');
    staid(ista) = scan{1};
    stacor(ista)  = scan{2};
    stacor_std(ista)  = scan{3};
    sta{ista}  = scan{4}{1};
end
fclose(fid);

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

