function [sta] = load_stations(STAPATH)

fid = fopen(STAPATH);
temp = fgetl(fid);
ii = 0;
while ~feof(fid)
    ii = ii+1;
    temp = fgetl(fid);
    vals = textscan(temp,'%s %f %f %f');
    sta.staname{ii} = vals{1}{1};
    sta.lats(ii) = vals{2};
    sta.lons(ii) = vals{3};
    sta.elevs(ii) = vals{4};
end
fclose(fid);

end

