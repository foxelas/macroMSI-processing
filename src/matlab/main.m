
%% main
close all;
%clc;

%Modify userSettings.csv for options
userSettingsFile = '..\..\conf\userSettings.csv'; 
setOpt(userSettingsFile);
readData; %% redo intial bg removal with labels  etc the images are corrupt

%actionGetMaps; 
