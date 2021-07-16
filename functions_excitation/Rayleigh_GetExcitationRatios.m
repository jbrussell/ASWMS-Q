
function [event] = Rayleigh_GetExcitationRatios(event,Periodlist,CARD,Path2Phv,Path2Eig,MaxN,MinorOrMajor)
%% Loops to Calculate Excitation ratios
%% for a specific period, and for a set of events

Depthlist = event.evdp;
Mrrlist = event.moment.m_rr;
Mttlist = event.moment.m_tt;
Mpplist = event.moment.m_pp;
Mrtlist = event.moment.m_rt;
Mrplist = event.moment.m_rp;
Mtplist = event.moment.m_tp;

% Loop over stations
for ista = 1:length(event.stadata)
    [~, Azimuthlist] = distance(event.evla,event.evlo,...
                                event.stadata(ista).stla,event.stadata(ista).stlo...
                                ,referenceEllipsoid('GRS80'));
    
    % Loop over overtones
    for currN = [0:1:MaxN]

        qfile = dir([Path2Phv,'/',CARD,'.s*to*.q']);
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
            Eigfname=[Path2Eig,'/',CARD,'.S',num2str(currN),'.mat'];
            tmpinfo=load(Eigfname);
            eig = tmpinfo.eig;
            pers_eig = [eig(:).per_want];
            [~,Iper] = min(abs(pers_eig-period));
            U=eig(Iper).u;
            Uderiv=eig(Iper).up;
            V=eig(Iper).v;
            Vderiv=eig(Iper).vp;
            r=eig(1).r;    

            if currN > 0 & MinorOrMajor
                wvgrpdx=2;
            else
                wvgrpdx=1;
            end

            % Loop over events
            for evnum = 1:1:length(Depthlist)

                 [ B_SourceAmp,B_SourcePhase,B_Complex_Rad_Pattern,B_Term1,...
                B_Term2,B_Term3 ] = ...
                GetRayleighSourceAmpandPhase([0:0.5:360],1000*Depthlist(evnum),period,...
                r, U, Uderiv,V,Vderiv,Mrrlist(evnum),Mttlist(evnum),...
                Mpplist(evnum),Mrtlist(evnum),Mrplist(evnum),...
                Mtplist(evnum),CurrC,wvgrpdx );  

                tmpmax = max(B_SourceAmp); maxB_SourceAmp = tmpmax(1);
                PeriodStruc(periodcounter).period(evnum) = period;
                PeriodStruc(periodcounter).mode(evnum,currN+1) = currN;
                PeriodStruc(periodcounter).Amp_max(evnum,currN+1) = maxB_SourceAmp;
                

                [ B_SourceAmp,B_SourcePhase,B_Complex_Rad_Pattern,B_Term1,...
                B_Term2,B_Term3 ] = ...
                GetRayleighSourceAmpandPhase(Azimuthlist(evnum),1000*Depthlist(evnum),period,...
                r, U, Uderiv,V,Vderiv,Mrrlist(evnum),Mttlist(evnum),...
                Mpplist(evnum),Mrtlist(evnum),Mrplist(evnum),...
                Mtplist(evnum),CurrC,wvgrpdx );  

                PeriodStruc(periodcounter).Amp(evnum,currN+1) = B_SourceAmp;
                PeriodStruc(periodcounter).Phase(evnum,currN+1) = B_SourcePhase;
                
                PeriodStruc(periodcounter).ratio_AmpMax(evnum,currN+1) = B_SourceAmp./ maxB_SourceAmp;

            end
        end  

    end
    % Save excitation structure
    event.stadata(ista).excitation = PeriodStruc;

    periodcounter=0;
    
    for period = Periodlist
        periodcounter=periodcounter+1;

        RawExcitation_Mat = PeriodStruc(periodcounter).Amp;

        if MaxN > 0
            for overtone_num = [1:1:MaxN]
                ExcitationRatio_Mat(:,overtone_num) = RawExcitation_Mat(:,overtone_num+1)./RawExcitation_Mat(:,1);
            end
        end

        event.stadata(ista).excitation(periodcounter).ratio_Amp1_0 = ExcitationRatio_Mat;
    end
end



