%% main 
close all; 
%clc; 

showImages = false;
saveImages = false; %true;
dataset = 'saitama_v8_min_region_bright';
%dataset = 'saitama_v8_min_region_dark';

tryReadData = false; %true;
    
%Set-up of options for running
options =  setOpt([], dataset, showImages, saveImages, tryReadData);   
readData; %% redo intial bg removal with labels  etc the images are corrupt 
%plotMeasuredSpectra(ID, Spectra, 1, options.saveOptions);
%load('D:\temp\Google Drive\titech\research\output\saitama_v8_min_region_bright\RecostructionComparison\reconstructionComparison.mat')
%ReflectanceEstimationParameterComparison;
%actionReflectanceEstimationComparison;
%actionLBP;
%actionRecostructSRGB;
%visualTool;