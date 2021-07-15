clear all;

x_max = -120; %-147; %-140; %-120; %-147; % degrees
x_min = -250; %-157; %-150; %-157; %-250; %-170; % degrees
lat = 16.600;
dx = 0.5; %0.5;
elev = 0;

lon = (x_min:dx:x_max);
lon = flip(lon);
I_lon = lon < -180;
lon(I_lon) = lon(I_lon) + 360;
num_sta = length(lon);

fid = fopen('synthetic_LRT.stn','w');
fprintf(fid,'%d,\n',num_sta);
for i = 1:num_sta
    fprintf(fid,'B%02d %10.3f %10.3f %6.3f\n',i,lat,lon(i),elev);
end
fclose(fid);