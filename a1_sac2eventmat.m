% This script is used to convert event based sac files into event mat which is the input for the following 
% scripts
%
%
clear;

setup_parameters;

comp = parameters.component;
isoverwrite = 1;
% dbpath = './sacdata/';
% eventfile = 'eventlist';
dbpath = parameters.dbpath;
eventfile = parameters.eventfile;
% outpath = './eventmat/';
workingdir = parameters.workingdir;
outpath = [workingdir,'eventmat/'];
minMw = parameters.minMw;
maxdepth = parameters.maxdepth;

if ~exist(outpath)
	mkdir(outpath);
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
		% build up event.stadata structure
		event.stadata(isac).stla = sac.STLA;
		event.stadata(isac).stlo = sac.STLO;
		event.stadata(isac).stel = sac.STEL;
		event.stadata(isac).dist = vdist(sac.STLA,sac.STLO,sac.EVLA,sac.EVLO)/1e3;
		event.stadata(isac).otime = otime+sac.B;
		event.stadata(isac).delta = sac.DELTA;
		event.stadata(isac).data = sac.DATA1;
		event.stadata(isac).cmp = sac.KCMPNM;
		event.stadata(isac).stnm = sac.KSTNM;
		event.stadata(isac).filename = sacfilename;
%		if sac.DELTA == 1
%			event.stadata(isac).cmp = 'LHZ';
%		end
    end
    if event.evdp<=maxdepth && event.Mw>=minMw
        matfilename = [outpath,char(eventids(ie)),'_',comp,'.mat'];
        save(matfilename,'event')
        disp(['Save to ',matfilename]);
    end
end % end of loop ie

