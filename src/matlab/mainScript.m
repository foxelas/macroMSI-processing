%% main 
close all; 
%clc; 
onOffOptions = {'OFF', 'ON'};

showImages = setAndNotify('show images', true);
saveImages = setAndNotify('save images', true);
tryReadData = setAndNotify('try read data', false);
dataset = setAndNotify('dataset', 'saitama_v9_bright_3class');

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