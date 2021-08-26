
function [eventcs] = Love_GetExcitationRatios(eventcs,Periodlist,CARD,Path2Phv,Path2Eig,MaxN,MinorOrMajor)
%% Loops to Calculate Excitation ratios
%% for a specific period, and for a set of events

Depthlist = eventcs.evdp;
Mrrlist = eventcs.moment.m_rr;
Mttlist = eventcs.moment.m_tt;
Mpplist = eventcs.moment.m_pp;
Mrtlist = eventcs.moment.m_rt;
Mrplist = eventcs.moment.m_rp;
Mtplist = eventcs.moment.m_tp;

% Loop over stations
for ista = 1:length(eventcs.autocor)
    [~, Azimuthlist] = distance(eventcs.evla,eventcs.evlo,...
                                eventcs.stlas(ista),eventcs.stlos(ista)...
                                ,referenceEllipsoid('GRS80'));
    
    % Loop over overtones
    for currN = [0:1:MaxN]

        qfile = dir([Path2Phv,'/',CARD,'.t*to*.q']);
        qfile = [Path2Phv,'/',qfile(1).name];
        mode = readMINEOS_qfile(qfile,currN);
        Tlist=mode.T;
        Clist=mode.phv;

        periodcounter = 0;

        % Loop over periods
        for period = Periodlist

            Tdiff = abs(Tlist-period);
            [mindiff,cdx]=min(Tdiff);
            CurrC = Clist(cdx);
            CurrC = deg2rad(km2deg(CurrC));

            periodcounter = periodcounter +1;
            % Load eigenfunction
            Eigfname=[Path2Eig,'/',CARD,'.T',num2str(currN),'.mat'];
            tmpinfo=load(Eigfname);
            eig = tmpinfo.eig;
            pers_eig = [eig(:).per_want];
            [~,Iper] = min(abs(pers_eig-period));
            W=eig(Iper).wl;
            Wderiv=eig(Iper).wp;
            r=eig(1).r;    

            if currN > 0 & MinorOrMajor
                wvgrpdx=2;
            else
                wvgrpdx=1;
            end

            % Loop over events
            for evnum = 1:1:length(Depthlist)

                 [ B_SourceAmp,B_SourcePhase,B_Complex_Rad_Pattern,B_Term1,...
                B_Term2 ] = ...
                GetLoveSourceAmpandPhase([0:0.5:360],1000*Depthlist(evnum),period,...
                r, W, Wderiv,Mrrlist(evnum),Mttlist(evnum),...
                Mpplist(evnum),Mrtlist(evnum),Mrplist(evnum),...
                Mtplist(evnum),CurrC,wvgrpdx );  

                tmpmax = max(B_SourceAmp); maxB_SourceAmp = tmpmax(1);
                PeriodStruc(periodcounter).period(evnum) = period;
                PeriodStruc(periodcounter).mode(evnum,currN+1) = currN;
                PeriodStruc(periodcounter).Amp_max(evnum,currN+1) = maxB_SourceAmp;
                

                [ B_SourceAmp,B_SourcePhase,B_Complex_Rad_Pattern,B_Term1,...
                B_Term2 ] = ...
                GetLoveSourceAmpandPhase(Azimuthlist(evnum),1000*Depthlist(evnum),period,...
                r, W, Wderiv,Mrrlist(evnum),Mttlist(evnum),...
                Mpplist(evnum),Mrtlist(evnum),Mrplist(evnum),...
                Mtplist(evnum),CurrC,wvgrpdx );  

                PeriodStruc(periodcounter).Amp(evnum,currN+1) = B_SourceAmp;
                PeriodStruc(periodcounter).Phase(evnum,currN+1) = B_SourcePhase;
                
                PeriodStruc(periodcounter).ratio_AmpMax(evnum,currN+1) = B_SourceAmp./ maxB_SourceAmp;

            end
        end  

    end
    % Save excitation structure
    eventcs.source(ista).excitation = PeriodStruc;

    periodcounter=0;
    
    for period = Periodlist
        periodcounter=periodcounter+1;

        RawExcitation_Mat = PeriodStruc(periodcounter).Amp;

        if MaxN > 0
            for overtone_num = [1:1:MaxN]
                ExcitationRatio_Mat(:,overtone_num) = RawExcitation_Mat(:,overtone_num+1)./RawExcitation_Mat(:,1);
            end
        end

        eventcs.source(ista).excitation(periodcounter).ratio_Amp1_0 = ExcitationRatio_Mat;
    end
end



