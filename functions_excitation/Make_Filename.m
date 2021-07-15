%%%%%%%%%%%%%%%%%%%%%%%
%%%% Make Output Filenames %%%%
%%%%%%%%%%%%%%%%%%%%%%%
if RayleighOrLove == 1
    Wave = 'R';
elseif RayleighOrLove == 0
    Wave = 'L';
end

if MinorOrMajor == 1
    Arc = 'Maj';
elseif MinorOrMajor == 0
    Arc = 'Mnr';
end


if ~exist('outdir','var')
    RawExcitationFname = [Wave '_' Arc '_Exct_T_' ...
        num2str(period) 's_output' Inputfilename(1:end-4)];
    ExcitationRatioFname = [Wave '_' Arc '_ExctRtio_T_'  ...
        num2str(period) 's_output' Inputfilename(1:end-4)];
else
    RawExcitationFname = [outdir '/' Wave '_' Arc '_Exct_T_' ...
        num2str(period) 's_output' Outputfilename];
    ExcitationRatioFname = [outdir '/' Wave '_' Arc '_ExctRtio_T_'  ...
        num2str(period) 's_output' Outputfilename];
end
