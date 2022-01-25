clear;
addpath('..')

setup_parameters_tpw;
comp = parameters.component;
periods = parameters.periods;
workingdir_tpw = parameters_tpw.workingdir;

phv_avg = [];
for ip = 1:length(periods)
    event = [];
    period = periods(ip);    
    control_file = [workingdir_tpw,'/','all.',num2str(round(period),'%03d')];
    phamp_file = [workingdir_tpw,'/','phampcor.',num2str(round(period),'%03d'),'.inp'];
    cfp = fopen(control_file,'r');
    stemp = fgetl(cfp);
    eventnum = sscanf(stemp,'%d');
    stanum=[];
    for ie = 1:eventnum
        stemp = fgetl(cfp);
        dtemp = sscanf(stemp,'%d %d');
        stanum(ie) = dtemp(1);
        for ista = 1:stanum(ie)
            stemp = fgetl(cfp);
        end
    end	
    fclose(cfp);

    dfp = fopen(phamp_file,'r');
    for ie = 1:eventnum
        stemp = fgetl(dfp);
        eventid = sscanf(stemp,'%d');
        if eventid ~=ie
            disp('something wrong!');
        end
        event(ie).id = eventid;
        for ista = 1:stanum(ie)
            stemp = fgetl(dfp);
            dtemp = sscanf(stemp,'%f %f');
            event(ie).sta(ista).bgtime = dtemp(1);
            event(ie).evid = num2str(dtemp(2));
            stemp = fgetl(dfp);
            dtemp = sscanf(stemp,'%f %f %f %f %f %f');
            event(ie).sta(ista).dist = dtemp(1);
            event(ie).sta(ista).azi = dtemp(2);
            event(ie).sta(ista).baz = dtemp(3);
            event(ie).sta(ista).deg = dtemp(4);
            event(ie).sta(ista).lat = dtemp(5);
            event(ie).sta(ista).lon = dtemp(6);
            stemp = fgetl(dfp);
            dtemp = sscanf(stemp,'%f %f');
            event(ie).sta(ista).amp = dtemp(1);
            event(ie).sta(ista).ph = dtemp(2);
        end
    end

    % normalize the phase

    alldist = [];
    allphs = [];
    allamps = [];
    for ie = 1:eventnum
        dists = [event(ie).sta.dist];
        phs = [event(ie).sta.ph];
        [mindist mindistid] = min(dists);
        dphs = phs - phs(mindistid);
        amps = [event(ie).sta.amp];
        damps = amps./amps(mindistid);
        ddists = dists - dists(mindistid);
        alldist = [alldist;ddists(:)];
        allphs = [allphs;dphs(:)];
        allamps = [allamps;damps(:)];
    end
    
    % Fit line to phase data to get average phase velocity    
    c = polyfit(alldist,allphs,1);
    phs_fit = polyval(c,alldist);
    phv_fit = 1./(c(1)*period);
    phv_avg(ip) = phv_fit;

    figure(23)
%     if ip == 1
%         clf;
%     end
    subplot(4,4,ip);
    sgtitle('Phase');
    hold on
    plot(alldist,allphs,'x');
    xdist = [min(alldist) max(alldist)];
    yphs = polyval(c,xdist);
    plot(xdist,yphs,'-r','linewidth',1.5);
    title([num2str(period),' s  ',num2str(phv_fit),' km/s']);
    % [x y] = ginput(2);
    % plot(x,y,'r');

    % phv = abs(diff(x))/abs(diff(y)*period)
    
    % Binned amplitudes
    dkm = 50;
    bins = min(alldist):dkm:max(alldist);
    cbins = bins(1:end-1) + dkm/2;
    amps_bin = [];
    for ibin = 1:length(bins)-1
        I = alldist>=bins(ibin) & alldist<bins(ibin+1);
        amps_bin(ibin) = nanmean(allamps(I));
    end
    
    figure(24)
%     if ip == 1
%         clf;
%     end
    subplot(4,4,ip);
    sgtitle('Amplitude');
    hold on
    plot(alldist,allamps,'x');
    plot(cbins,amps_bin,'-or');
    plot([min(alldist),max(alldist)],[1 1],'--k');
    title([num2str(period),' s']);
end

%%
figure(25); hold on;
plot(periods,phv_avg,'-ob');
xlabel('Period (s)');
ylabel('Phase Velocity (km/s)');

path2qfile = '../qfiles/pa5_5km.s0to66.q';
if exist(path2qfile,'file')==2
    mineos = readMINEOS_qfile(path2qfile,0);
    phv_mineos = interp1(mineos.T,mineos.phv,periods);
    plot(periods,phv_mineos,'-r');
end

