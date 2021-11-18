%% This script is used for convert the output of matGSDF into input files for Two-Plane-Wave fortran
%  code. 
%  Also, generate stanumlist (list of station inputs)
%
% This version uses source excitation estimates from MINEOS to avoid nodes in 
% the radiation pattern. Must first run MINEOS to generate eigenfunction files
% - ../MINEOS/run_MINEOS/a1_run_mineos_check.m
% - ../MINEOS/run_MINEOS/a2_mk_kernels.m
% - ../MINEOS/run_MINEOS/a4_pull_eigenfuncs.m (Fund. and 1st mode)
% Then attach excitation values to CSmeasure by running ../a5_a_attach_excitation.m
% and copy new CSmeasure files into the TPW directory.
%
%  JBR 10/17/2021: modified for TPW version that includes anisotropy and
%  attenuation
%
clear;

setup_parameters_tpw;
comp = parameters.component;
periods = parameters.periods;
refvs = parameters_tpw.refphv; %[3.53 3.61 3.7 3.76 3.82 3.85];
stalist = parameters_tpw.stalist;
min_excitation_ratio = parameters_tpw.min_excitation_ratio; % minimum source excitation ratio (A/A_max) to consider. This avoids nodes in radiation pattern

% input path
workingdir = [parameters.workingdir];
eventcs_path = [workingdir,'CSmeasure/'];
eikonal_output_path = [workingdir,'eikonal/'];
% eventcs_path = './CSmeasure/';
% eikonal_output_path = './eikonal/';

% Parameters
workingdir_tpw = parameters_tpw.workingdir;

min_sta_num = parameters_tpw.min_sta_num;
nfreq = parameters_tpw.nfreq;
iterlimit = parameters_tpw.iterlimit;
divfac = parameters_tpw.divfac;
divfac_azi = parameters_tpw.divfac_azi;
dampvel = parameters_tpw.dampvel;
dampaniso = parameters_tpw.dampaniso;
refalpha = parameters_tpw.refalpha;
dampalpha = parameters_tpw.dampalpha;

stnms_all = {};
for ip = 1:length(periods)
    % output files
    period = periods(ip);
    control_file = [workingdir_tpw,'/','all.',num2str(period,'%03d')];
    ampph_file = ['phampcor.',num2str(period,'%03d'),'.inp'];

    % input parameters
    freq = 1/periods(ip);
    detailoutput = ['detail.',num2str(period,'%03d'),'.inp'];
    summaroutput = ['summar.',num2str(period,'%03d'),'.inp'];
    gridnodes = ['gridinp.dat'];
    varfile = ['covar.',num2str(period,'%03d'),'.inp'];
    mavamp = ['mavamp.',num2str(period,'%03d'),'.inp'];
    fvelarea0 = ['temp.',num2str(period,'%03d'),'.inp'];
    velarea = ['velarea.',num2str(period,'%03d'),'.inp'];
    unifvel = refvs(ip);
    kernelfiledir = dir([workingdir_tpw,'/','sensspec',num2str(period,'%03d'),'s*']);
    kernelfile = [kernelfiledir(1).name];
    outgrd = ['outgrd.',num2str(period,'%03d'),'.inp'];
    outvel = ['outvel.',num2str(period,'%03d'),'.txt'];
    outazi = ['outazi.',num2str(period,'%03d'),'.txt'];
    outstacor = ['outstacor.',num2str(period,'%03d'),'.txt'];
    outalpha = ['outalpha.',num2str(period,'%03d'),'.txt'];

    % read in bad station list, if existed
    if exist('badsta.lst')
        badstnms = textread('badsta.lst','%s');
        disp('Found Bad stations:')
        disp(badstnms)
    end

    % gather necessary information from GSDF output
    eikonal_matfiles = dir([eikonal_output_path,'/*eikonal_',comp,'.mat']);
    goodeventnum = 0;
    for ie=1:length(eikonal_matfiles)
        temp = load([eikonal_output_path,eikonal_matfiles(ie).name]);
        eventphv = temp.eventphv;
        disp(eventphv(ip).id);
        stlas = eventphv.stlas;
        stlos = eventphv.stlos;
        stnms = eventphv.stnms;
        if exist('badstnms','var')
            badstaids = find(ismember(stnms,badstnms));
        else
            badstaids = [];
        end
        tp = eventphv(ip).traveltime;
        goodstaind = find(~isnan(tp));
        if ~isempty(badstaids)
            badind = find(ismember(goodstaind,badstaids));
            goodstaind(badind) = [];
        end
        goodstanum = length(goodstaind);
        if goodstanum > min_sta_num
            temp = load([eventcs_path,'/',eventphv(ip).id,'_cs_',comp,'.mat']);
            eventcs = temp.eventcs;
            
            amp_ratios = [];
            for i = 1:length(goodstaind)
                ista = goodstaind(i);
                amp_ratios(i) = eventcs.source(ista).excitation(ip).ratio_AmpMax(1);
            end
            if mean(amp_ratios) < min_excitation_ratio
                continue
            end
            
            goodeventnum = goodeventnum + 1;
            twpevents(goodeventnum).evla = eventphv(ip).evla;
            twpevents(goodeventnum).evlo = eventphv(ip).evlo;
            twpevents(goodeventnum).stlas = eventphv(ip).stlas(goodstaind);
            twpevents(goodeventnum).stlos = eventphv(ip).stlos(goodstaind);
            twpevents(goodeventnum).stlos(twpevents(goodeventnum).stlos<0) = twpevents(goodeventnum).stlos(twpevents(goodeventnum).stlos<0) + 360;
            twpevents(goodeventnum).stnms = eventphv(ip).stnms(goodstaind);
            twpevents(goodeventnum).tp = eventphv(ip).traveltime(goodstaind);
            twpevents(goodeventnum).id = eventphv(ip).id;
            for i = 1:length(goodstaind)
                ista = goodstaind(i);
                twpevents(goodeventnum).amp(i) = eventcs.autocor(ista).amp(ip);
            end
        end
    end
    
    % Save unique station names from good events
    stnms_all{ip} = unique([twpevents.stnms]);

    % Calculate useful informations for TPW output
    for ie=1:length(twpevents)
        stlas = twpevents(ie).stlas;
        stlos = twpevents(ie).stlos;
        evla = twpevents(ie).evla;
        evlo = twpevents(ie).evlo;
        degs = km2deg(distance(stlas,stlos,evla,evlo,referenceEllipsoid('GRS80'))/1000);
        twpevents(ie).degs = degs;
        twpevents(ie).dists = deg2km(degs);
        twpevents(ie).azi = azimuth(evla,evlo,stlas,stlos,referenceEllipsoid('GRS80'));
        twpevents(ie).baz = azimuth(stlas,stlos,evla,evlo,referenceEllipsoid('GRS80'));
    end

    % output the files
    cfp = fopen(control_file,'w');
    afp = fopen([workingdir_tpw,'/',ampph_file],'w');

    fprintf(cfp,'%d\n',length(twpevents));
    for ie=1:length(twpevents)
        fprintf(cfp,'%d   %d\n',length(twpevents(ie).stlas),ie);
        fprintf(afp,' %d\n',ie);
        maxtp = max([twpevents(ie).tp]);
        mintp = min([twpevents(ie).tp]);
        meanamp = mean([twpevents(ie).amp].^0.5);
        for ista=1:length(twpevents(ie).stlas)
            fprintf(cfp,'%s\n',char(twpevents(ie).stnms(ista)));
            fprintf(afp,' 0.0\n');
            dist = twpevents(ie).dists(ista);
            azi = twpevents(ie).azi(ista);
            baz = twpevents(ie).baz(ista);
            stla = twpevents(ie).stlas(ista);
            stlo = twpevents(ie).stlos(ista);
            deg = twpevents(ie).degs(ista);
            amp = twpevents(ie).amp(ista).^0.5;
            stnm = char(twpevents(ie).stnms(ista));
            amp_mean = amp./meanamp;
            tp = twpevents(ie).tp(ista);
    %		ph = (maxtp-tp)./periods(ip)/2/pi;
            ph = (tp-mintp)./periods(ip);
            fprintf(afp,' %f  %f  %f  %f  %f  %f\n',dist,azi,baz,deg,stla,stlo);
            fprintf(afp,' %d  %f\n',amp,ph);
        end
    end

    fprintf(cfp,'%d\n',nfreq); % nfreq - number of requencies (always 1)
    fprintf(cfp,'%f\n',freq); % freq - frequency of interest
    fprintf(cfp,'%s\n',detailoutput); % foutput - detail.???.inp
    fprintf(cfp,'%s\n',summaroutput); % fsummary - summar.???.inp
    fprintf(cfp,'%s\n',gridnodes); % finvrsnodes - gridinp
    fprintf(cfp,'%s\n',ampph_file); % fftinput - phampcor.???.inp 
    fprintf(cfp,'%s\n',varfile); % fvariance - covar.???.inp
    fprintf(cfp,'%s\n',mavamp); % fmaxavamp - mavamp.???.inp
    fprintf(cfp,'%s\n','dtemp'); % ftemp - ftemp.???.inp
    fprintf(cfp,'%s\n','ddd'); % fvelarea0 - temp.???.inp
    fprintf(cfp,'%s\n',stalist); % fstalist
    fprintf(cfp,'%f\n',unifvel); % unifvel - reference velocity
    fprintf(cfp,'%d %f %f %f %f\n',iterlimit,dampvel,dampaniso,divfac,divfac_azi); % interlimit, dampvel, dampaniso, divfac
    fprintf(cfp,'%s\n',kernelfile); % sensfn - sensitivity kernels
    fprintf(cfp,'%s\n','velout.dum'); % fvelout - output velocity file
    fprintf(cfp,'%s\n',outgrd); % startvel
    fprintf(cfp,'%s\n','resolve.dum'); % fresdiag - Resolution
    fprintf(cfp,'%s\n','gridphase.dum'); % fendvel
    fprintf(cfp,'%s\n',outvel); % output phase velocity
    fprintf(cfp,'%s\n',outazi); % output azimuthal anisotropy
    fprintf(cfp,'%s\n',outstacor); % output station correction terms
    fprintf(cfp,'%s\n',outalpha); % output attenuation
    fprintf(cfp,'%f\n',refalpha); % refgamma - reference alpha
    fprintf(cfp,'%f\n',dampalpha); % dampgamma - alpha damping

    fclose(cfp);
    fclose(afp);

end

% Build station file
stnms = unique([stnms_all{:}]);
stafile = [workingdir_tpw,'/','stanumlist'];
fid = fopen(stafile,'w');
for ista = 1:length(stnms)
    fprintf(fid,'%s\n',stnms{ista});
end
fprintf(fid,'%s\n','nope'); % not sure if this is necessary, but Zhitu's version included it
fclose(fid);