function [moment] = load_moment_tensor(fname)
% Load moment tensor values from CMT idagrn file
fid = fopen(fname,'r');
moment.evid = sscanf(fgetl(fid),'%s');
temp = sscanf(fgetl(fid),'%f %f %f');
moment.lat = temp(1);
moment.lon = temp(2);
moment.depth_km = temp(3);
moment.mult_fac = sscanf(fgetl(fid),'%f');
moment.m_rr = sscanf(fgetl(fid),'%f');
moment.m_tt = sscanf(fgetl(fid),'%f');
moment.m_pp = sscanf(fgetl(fid),'%f');
moment.m_rt = sscanf(fgetl(fid),'%f');
moment.m_rp = sscanf(fgetl(fid),'%f');
moment.m_tp = sscanf(fgetl(fid),'%f');    
fclose(fid);

end

