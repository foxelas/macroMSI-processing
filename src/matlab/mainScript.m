%% main 
close all; clc; 

showImages = true;
saveImages = false; %true;
dataset = 'saitama_v8_min_region_bright';
tryReadData = false; %true;
    
%Set-up of options for running
options =  setOpt([], dataset, showImages, saveImages, tryReadData);   
readData; %% redo intial bg removal with labels  etc the images are corrupt 
options.action = 'Refest_Preset_plusrgb';
ReflectanceEstimationParameterComparison;
%actionReflectanceEstimationComparison;


