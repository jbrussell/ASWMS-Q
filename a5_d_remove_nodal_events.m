% Remove nodel events
% Remove events for which the array region is near a node in the radiation
% pattern. Excitation ratio is defined as A / A_max for a given
% station. A station that is near a node will have a low excitation ratio.
% This script removes a measurement if either of the pair of stations A/A_max is 
% less than the threshold set by A_Amax_minthresh.
%
% jbrussell 11/23/2021
%
clear

% setup parameters
setup_parameters
setup_ErrorCode

A_Amax_minthresh = 0.6; % Threshold for minimum excitation ratio

workingdir = parameters.workingdir;
receiverterms_path = [workingdir];
% input path
eventcs_path = [workingdir,'CSmeasure/'];
nodal_eventcs_path = [workingdir,'CSmeasure_nodal/'];
comp = parameters.component;
periods = parameters.periods;

csmatfiles = dir([eventcs_path,'/*cs_',comp,'.mat']);
ikill = 0;
for ie = 1:length(csmatfiles)
	% read in data and set up useful variables
	temp = load([eventcs_path,csmatfiles(ie).name]);
	eventcs =  temp.eventcs;
	disp(eventcs.id)
    
    all_ratios = [];
    ii = 0;
    for ip = 1:length(periods)
        for ista = 1:length(eventcs.autocor)
            ii = ii + 1;
            all_ratios(ip,ista) = eventcs.source(ista).excitation(ip).ratio_AmpMax(1);
        end
    end
    
    for imeas = 1:length(eventcs.CS)
        ista1 = eventcs.CS(imeas).sta1;
        ista2 = eventcs.CS(imeas).sta2;
        Ibad_pers = find(all_ratios(:,ista1)<A_Amax_minthresh | all_ratios(:,ista2)<A_Amax_minthresh);
        eventcs.CS(imeas).isgood(Ibad_pers) = ErrorCode.near_node;
    end
    
    for ista = 1:length(eventcs.autocor)
        Ibad_pers = find(all_ratios(:,ista)<A_Amax_minthresh);
        eventcs.autocor(ista).exitflag(Ibad_pers) = ErrorCode.near_node;
    end
    
    save([eventcs_path,csmatfiles(ie).name],'eventcs');

end % end of loop ie