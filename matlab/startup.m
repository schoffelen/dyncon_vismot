function startup

if ispc
    addpath('P:\3011085.03\scripts\fieldtrip');
    addpath('P:\3011085.03\scripts\fieldtrip\qsub');
    addpath('P:\3011085.03\scripts\fieldtrip\external\brewermap');
    addpath('M:\MATLAB');
    addpath('M:\MATLAB\cellfunction');
    cd('P:\3011085.03\scripts\dyncon_vismot\matlab');
else
    addpath('/project/3011085.03/scripts/fieldtrip');
    addpath('/project/3011085.03/scripts/fieldtrip/qsub');
    addpath('/project/3011085.03/scripts/fieldtrip/external/brewermap');
    
    % project folder
    addpath('/project/3011085.03/scripts/dyncon_vismot/matlab/');
    addpath('/project/3011085.03/scripts/dyncon_vismot/matlab_old/');
    
    % cell function
    addpath('/project/3011085.03/scripts/cellfunction'); 
    addpath('/project/3011085.03/scripts/fieldtrip/external/brewermap');
    
    % Canlab Modulation toolbox
    addpath(genpath('/project/3011085.03/scripts/ModulationToolbox/'));
    addpath(genpath('/project/3011085.03/scripts/CanlabCore/'));
    
    % alternative boxplot
    addpath /project/3011085.03/scripts/IoSR-Surrey-MatlabToolbox-4bff1bb/
        
    cd('/project/3011085.03/scripts/dyncon_vismot/matlab/');
end
ft_defaults

end