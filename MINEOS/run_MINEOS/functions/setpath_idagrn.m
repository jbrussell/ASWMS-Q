% Set path for idagrn executables
% 
% JBR, 2018

PATH = getenv('PATH');
setenv('SACHOME','/opt/local/sac');
setenv('PATH',[PATH,':/opt/local/sac/bin/']);
setenv('SACAUX','/opt/local/sac/aux');
