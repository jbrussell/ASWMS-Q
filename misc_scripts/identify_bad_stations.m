% Script to identify bad stations by counting the number of events each station has "good" recordings for. Bad stations are those that have good recordings for less than 25% of events.

clear; close all;

proj = 'OUT';

path2CSmeasure = ['./',proj,'/CSmeasure/'];

files = dir(path2CSmeasure);
files = files(3:end);

stnms = {};
sta_goodcount = [];
for iev = 1:length(files)
    filename = files(iev).name;
    temp = load([path2CSmeasure,'/',filename]);
    eventcs = temp.eventcs; clear temp;
    
    for ista = 1:length(eventcs.stnms)
        if isempty(find(strcmp(stnms,eventcs.stnms{ista})))
            stnms{end+1} = eventcs.stnms{ista};
            sta_goodcount(end+1) = 1;
        else
            Ista = strcmp(stnms,eventcs.stnms{ista});
            sta_goodcount(Ista) = sta_goodcount(Ista) + 1;
        end
    end
end

[stnms,isrt] = sort(stnms);
sta_goodcount = sta_goodcount(isrt);

%% Plot number of recordings per station
figure(1);
bh = bar(categorical(stnms),sta_goodcount);
title('Events per station: red = stations with good recordings for < 25% of events');
ylabel('Number of good recordings');
% Color bad stations red
Ibad = find(sta_goodcount < length(files)*0.25); % stations that see less than 25% of events
bh.FaceColor = 'flat'; % allow per-bar coloring
bh.CData(Ibad,:) = repmat([1 0 0],length(Ibad),1);