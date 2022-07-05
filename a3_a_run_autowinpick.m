% Scripts to run the auto_win_pick function for all the events and generate a old version "events" file
clear;
%plot native;

isfigure = 1;
is_overwrite = 1;

% Setup parameters
setup_parameters

eventmatpath = './eventmat/';
outwinpath = './winpara/';
% workingdir = parameters.workingdir;
% eventmatpath = [workingdir,'eventmat/'];
% outwinpath = [workingdir,'winpara/'];

if ~exist(outwinpath,'dir')
	mkdir(outwinpath)
end

% Setup Error Codes for Bad data
setup_ErrorCode

periods = parameters.periods;
comp = parameters.component;

matfiles = dir([eventmatpath,'/*_',comp,'.mat']);
for ie = 1:length(matfiles)

	clear event eventcs CS
	% read in the events information
	temp = load([eventmatpath,matfiles(ie).name]);
	event = temp.event;
	disp(event.id)

	if ~is_overwrite
		filename = [outwinpath,'/',event.id,'_',comp,'.bad'];
		if exist(filename,'file')
			disp(['Exist: ',filename,' Skip!']);
			continue;
		end
		filename = [outwinpath,'/',event.id,'_',comp,'.win'];
		if exist(filename,'file')
			disp(['Exist: ',filename,' Skip!']);
			continue;
		end
	end

	% set up some useful arrays
    if ~isfield(event,'stadata')
        continue;
    end
    stlas = [event.stadata(:).stla];
    stlos = [event.stadata(:).stlo];
    stnms = {event.stadata(:).stnm};
    dists = [event.stadata(:).dist];
    

	% automatically select the signal window by using ftan method
	disp('Start to picking the window');
	tic
		[winpara event] = auto_win_select(event);
	toc
	if isfigure && length(find([event.stadata.isgood]>0))>1
		plot_win_select(event,periods,winpara);
	end
	if length(winpara) ~= 4
		filename = [outwinpath,'/',event.id,'_',comp,'.bad'];
		fp = fopen(filename,'w');
		fprintf(fp,'%f\n',0);
		fclose(fp);
		continue;
	end
	filename = [outwinpath,'/',event.id,'_',comp,'.win'];
	fp = fopen(filename,'w');
	fprintf(fp,'%s %f %f %f %f\n',event.id,winpara(1),winpara(2),winpara(3),winpara(4));
	fclose(fp);
	event.winpara = winpara;
	save([eventmatpath,matfiles(ie).name],'event');
	disp(winpara);
	
end % end of event
