% ************************************************************************
% This toolbox is a temporary package from the work of Louis Korczowski @GIPSA-lab
% In this state, it is only made from personal testing and not public release.
%
% *** History: 2016-05-04
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2016
% *** Related work:  M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis", 2016.
%
% see also: ACSTP, script_ACSTP_test


p=mfilename('fullpath'); %find the path from this script
p=p(1:end-10);%remove 'installer'

% install the external
DIR=dir([p '/external/']);
DIR=DIR(3:end);

for indE=1:length(DIR)
    if DIR(indE).isdir==1
        pathEx=[p '/external/'];
        pathExDir=[pathEx DIR(indE).name];
        addpath(pathExDir);
        DIR2=dir(pathExDir);DIR2=DIR2(3:end);
        for indE2=1:length(DIR2)
            if strcmp(DIR2(indE2).name,'installer.m')
                run( [pathExDir '/' (DIR2(indE2).name)])
                p=mfilename('fullpath'); %find the path from this script
                p=p(1:end-10);%remove 'installer'
            end
        end
    end
    
end
addpath(p);
addpath( genpath([p '/data/']) );

addpath( genpath([p '/lib/']) );
addpath( genpath([p '/scripts/']));

disp('LK TOOLBOX successfully activated')
help installer
clear p

