% Manually fix automatic windows
clear;
%plot native

% eventmat_files = dir('eventmat/*.mat');
setup_parameters
workingdir = parameters.workingdir;
eventmat_files = dir([workingdir,'eventmat/*.mat']);

for ie=1:length(eventmat_files)
% 	load(fullfile('eventmat',eventmat_files(ie).name));
    load(fullfile(workingdir,'eventmat',eventmat_files(ie).name));
	if isfield(event,'winpara')
		if ~isfield(event,'isgood')
			event.isgood = 1;
		end
		disp(eventmat_files(ie).name);
        disp('p for pick window, g for isgood, q for next event');
        disp('c for cheat sheet');
        event = exam_winpara_phase(event);
	end
% 	save(fullfile('eventmat',eventmat_files(ie).name),'event');
    save(fullfile(workingdir,'eventmat',eventmat_files(ie).name),'event');
end

