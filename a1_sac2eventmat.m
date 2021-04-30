% This script is used to convert event based sac files into event mat which is the input for the following 
% scripts
%
%
clear;

setup_parameters;

comp = parameters.component;
isoverwrite = 1;
iscleanevs = 1;
isplotsnr = 0;
% dbpath = './sacdata/';
% eventfile = 'eventlist';
dbpath = parameters.dbpath;
eventfile = parameters.eventfile;
% outpath = './eventmat/';
workingdir = parameters.workingdir;
outpath = [workingdir,'eventmat/'];
minMw = parameters.minMw;
maxdepth = parameters.maxdepth;
max_dist_tol = parameters.max_dist_tol;
snr_tol = parameters.snr_tol;

if ~exist(outpath)
	mkdir(outpath);
end
if iscleanevs
    system(['rm ',outpath,'*']);
end

eventids = textread([dbpath,eventfile],'%s');
for ie = 1:length(eventids)
	matfilename = [outpath,char(eventids(ie)),'_',comp,'.mat'];
	if ~isoverwrite && exist(matfilename)
		disp(['Exist ',matfilename,', Skip!']);
		continue;
	end
	clear event
	datapath = [dbpath, char(eventids(ie)),'/'];
	disp(datapath);
	saclist = dir([datapath,'*',comp,'.sac']);
    if isempty(saclist)
        continue
    end
    iisac = 0;
    is_skip_mag_dep = 0; % Initialize to 0
	for isac = 1:length(saclist)
		% read sac file
		sacfilename = [datapath,saclist(isac).name];
		sac = readsac(sacfilename);
		days = (datenum(sac.NZYEAR-1,12,31)+sac.NZJDAY-datenum(1970,1,1));
		otime = days*24*3600 + sac.NZHOUR*3600+sac.NZMIN*60+sac.NZSEC+sac.NZMSEC/1000;
		% initial the event information by the first sac file
		if isac == 1
            if sac.EVDP/1000 > maxdepth || sac.MAG < minMw
                is_skip_mag_dep = 1;
                disp(['Magnitude or depth exceeded for ',matfilename,', Skip!']);
                break;
            end
			event.evla = sac.EVLA;
			event.evlo = sac.EVLO;
            event.evdp = sac.EVDP/1000;
			event.otime = otime;
			event.dbpath = datapath;
			event.id = char(eventids(ie));
			event.otimestr = datestr(datenum(sac.NZYEAR-1,12,31,sac.NZHOUR,sac.NZMIN,sac.NZSEC)+sac.NZJDAY);
            event.Mw = sac.MAG;
        end
		
		% wbh check for timing errors
        %disp(abs(event.otime - otime))
        if abs(event.otime - otime) > 1
            error(['Data start time more than 1 s different from event time for ',matfilename]);
        end
		
        % calculate snr
        snr = calc_SNR(sac,parameters,isplotsnr);
        if snr <= snr_tol || isnan(snr)
%             disp('low snr... skipping')
            continue
        end
        iisac = iisac+1;
		% build up event.stadata structure
		event.stadata(iisac).stla = sac.STLA;
		event.stadata(iisac).stlo = sac.STLO;
		event.stadata(iisac).stel = sac.STEL;
		event.stadata(iisac).dist = vdist(sac.STLA,sac.STLO,sac.EVLA,sac.EVLO)/1e3;
		event.stadata(iisac).otime = otime+sac.B;
		event.stadata(iisac).delta = sac.DELTA;
		event.stadata(iisac).data = sac.DATA1;
		event.stadata(iisac).cmp = sac.KCMPNM;
		event.stadata(iisac).stnm = sac.KSTNM;
		event.stadata(iisac).filename = sacfilename;
        event.stadata(iisac).snr = snr;
%		if sac.DELTA == 1
%			event.stadata(isac).cmp = 'LHZ';
%		end
    end
    try 
        if mean([event.stadata(:).dist]) > max_dist_tol
            disp(['Max distance exceeded for ',char(eventids(ie)),', Skip!']);
            continue;
        end
    end
    if ~is_skip_mag_dep
        matfilename = [outpath,char(eventids(ie)),'_',comp,'.mat'];
        save(matfilename,'event')
        disp(['Save to ',matfilename]);
    end
end % end of loop ie

