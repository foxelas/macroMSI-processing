%% main 
close all; clc; 

showImages = true;
saveImages = true;
dataset = 'saitama_v8_min_region_bright';
tryReadData = true;
    
%Set-up of options for running
options =  setOpt([], dataset, showImages, saveImages, tryReadData);   
%readData;
options.action = 'Refest_Preset_plusrgb';
actionReflectanceEstimationComparison;


