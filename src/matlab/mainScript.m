%% main 
close all; 
%clc; 

showImages = setAndNotify('show images', true);
saveImages = setAndNotify('save images', true);
tryReadData = setAndNotify('try read data', false);
dataset = setAndNotify('dataset', 'saitama_v9_bright_3class'); 

options =  setOpt([], dataset, showImages, saveImages, tryReadData);   
readData; %% redo intial bg removal with labels  etc the images are corrupt 
% visualizePOIs;
% actionSOM;

[WMeasured, score, latent, explained] = dimensionReduction('PCA', Gun, double(labelsun));

% %plotMeasuredSpectra(ID, Spectra, 1, options.saveOptions);
% %ReflectanceEstimationParameterComparison;
% %actionReflectanceEstimationComparison;
% %actionLBP;
% %actionRecostructSRGB;
% %visualTool;