% Remove radiation pattern
% Radiation pattern is calculated and saved in a3_attach_excitation.m
%
% A_corr = A / A_radiation
clear

% setup parameters
setup_parameters

workingdir = parameters.workingdir;
receiverterms_path = [workingdir];
% input path
eventcs_path = [workingdir,'CSmeasure/'];
comp = parameters.component;
periods = parameters.periods;

csmatfiles = dir([eventcs_path,'/*cs_',comp,'.mat']);
for ie = 1:length(csmatfiles)
%for ie = 30
	clear ampgrad 
	% read in data and set up useful variables
	temp = load([eventcs_path,csmatfiles(ie).name]);
	eventcs =  temp.eventcs;
	disp(eventcs.id)
	
	% Loop through the periods
	for ip = 1:length(periods)
        
        % Divide amplitude data by |sind(delta)|^-0.5 to remove geometrical spreading
		for ista = 1:length(eventcs.autocor)
            
            % Get radiation pattern info
            mode = 0; % 0 = fundamental mode; 1 = first overtone, etc...
            AmpFactor = eventcs.source(ista).excitation(ip).ratio_AmpMax(mode+1);
            
			if eventcs.autocor(ista).exitflag(ip)>0
				amps(ista) = eventcs.autocor(ista).amp(ip);
			else
				amps(ista) = NaN;
            end
            
            % A_corr^2 = A^2 / Amp_radiation^2
            amp_corr = eventcs.autocor(ista).amp(ip) / AmpFactor.^2;
            eventcs.autocor(ista).amp(ip) = amp_corr;            
        end
    end
    
    save([eventcs_path,csmatfiles(ie).name],'eventcs');

end % end of loop ie
