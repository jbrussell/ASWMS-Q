% Function to read idagrn source excitation file
% JBR 9/2021

function [excite] = load_excitation_asc(path2asc,MODE)

% Read .asc file
dat = {};
imode = MODE+1;
com = ['awk ''{ if ($1 ==',num2str(imode-1),') print $0}'' ',path2asc];
[log3, dat{imode}] = system(com);
dat{imode} = str2num(dat{imode});
nn =  dat{imode}(:,1); % mode number
ll =  dat{imode}(:,2); % angular order
w = dat{imode}(:,3); % angular frequency (rad/s)
amp =  dat{imode}(:,4); % ak amplitude kernel
phase = dat{imode}(:,5); % phi2 source phase
phip =  dat{imode}(:,6); % phip ? (zero)
ae =  dat{imode}(:,7); % ae ?

excite.n = nn;
excite.l = ll;
excite.w = w;
excite.per = 2*pi./w;
excite.amp = amp;
excite.phase = phase;
excite.phip = phip;
excite.ae = ae;

end