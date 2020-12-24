
%% main
close all;
clear all;
clc;

%Modify userSettings.csv for options
userSettingsFile = '..\..\conf\userSettings.csv';
setOpt(userSettingsFile);

%% To keep settings in the workspace
configuration = load('configuration.mat');
readData;

%actionGetMaps;
