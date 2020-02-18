
%% main
close all;
%clc;

showImages = setAndNotify('show images', true);
saveImages = setAndNotify('save images', false);
tryReadData = setAndNotify('try read data', false);
dataset = setAndNotify('dataset', 'saitama_v9_bright_3class');

options = setOpt([], dataset, showImages, saveImages, tryReadData);
readData; %% redo intial bg removal with labels  etc the images are corrupt
% visualizePOIs;
% actionSOM;
% dimredscript;
% actionReadH5;

% %plotMeasuredSpectra(ID, Spectra, 1, options.saveOptions);
% %ReflectanceEstimationParameterComparison;
% %actionReflectanceEstimationComparison;
% %actionLBP;
% %actionRecostructSRGB;
% %visualTool;
createOpticalDensityMaps; 