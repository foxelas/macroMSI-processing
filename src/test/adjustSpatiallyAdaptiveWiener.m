close all; clc;
dataset = 'saitama_v2';
action = lower('ReflectanceEstimationSimple');
skipLoading = false;
showImages = false;
saveImages = false;

saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
options = struct('tryReadData', false, 'dataset', dataset, 'action', action, ...
'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
'showImages', showImages, 'saveOptions', saveOptions);
setup;
fminsearchOptions = optimset('MaxIter', 10, 'PlotFcns',@optimplotfval);
fileID = fopen('..\logs\optimization.txt', 'a');

tic;
%% Optimize noise model "givenSNR"
[minSNR,minVal] = fminsearch(@estimationGivenSNR, 10, fminsearchOptions);
outputLog = sprintf('Given SNR:: MinRMSE=%.5f, SNR=%d\n', minVal, minSNR)
fprintf(fileID, outputLog);

%% Optimize noise model "White Gaussian"
[minNoiseOrder,minVal] = fminsearch(@estimationWhiteGaussian, 2, fminsearchOptions);
outputLog = sprintf('White gaussian:: MinRMSE=%.5f, NoiseOrder=10^(-%d)\n', minVal, minNoiseOrder)
fprintf(fileID, outputLog);

%% Optimize noise model "Independent noise"
[minNoiseOrder,minVal] = fminsearch(@estimationIndependentNoise, 2, fminsearchOptions);
outputLog = sprintf('Independent noise:: MinRMSE=%.5f, NoiseOrder=10^(-%d)\n', minVal, minNoiseOrder)
fprintf(fileID, outputLog);
clear('rmse');

%% Optimize noise model "From Olympus"
options.noiseType = 'fromOlympus';
actionReflectanceEstimationComparison;
outputLog = sprintf('Simple Wiener:: MinRMSE=%.5f, Noise from Olympus\n', mean(rmse, 2))
fprintf(fileID, outputLog);
clear('rmse');

%% Optimize noise model "Spatial From Olympus"
options.noiseType = 'spatial';
actionReflectanceEstimationComparison;
outputLog = sprintf('Spatially adaptive Wiener:: MinRMSE=%.5f, Noise from Olympus\n', mean(rmse, 2))
fprintf(fileID, outputLog);
clear('rmse');

%% Optimize noise model "Spatial Indipendent"
estimationSpatiallyAdaptive
[sigma,minVal] = fminsearch(@estimationSpatiallyAdaptive, [10^(-4), 10^(-4)], fminsearchOptions);
outputLog = sprintf('Spatially adaptive Wiener:: MinRMSE=%.5f, Sigma1=%.4f, Sigma2=%.4f\n', minVal, sigma(1), sigma(2))
fprintf(fileID, outputLog);

fclose(fileID);
t = toc;
fprintf('Action elapsed %.2f mins.\n', t / 60);
clear variables; 

function rmseCur = estimationGivenSNR(x)
    dataset = 'saitama_v2';
    action = lower('ReflectanceEstimationSimple');
    skipLoading = false;
    showImages = false;
    saveImages = false;

    saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
    options = struct('tryReadData', false, 'dataset', dataset, 'action', action, ...
    'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
    'showImages', showImages, 'saveOptions', saveOptions);
    setup;
    options.noiseType = 'givenSNR';
    options.snr = x;
    actionReflectanceEstimationComparison;
    rmseCur = mean(rmse, 2);
end

function rmseCur = estimationWhiteGaussian(x)
    dataset = 'saitama_v2';
    action = lower('ReflectanceEstimationSimple');
    skipLoading = false;
    showImages = false;
    saveImages = false;

    saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
    options = struct('tryReadData', false, 'dataset', dataset, 'action', action, ...
    'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
    'showImages', showImages, 'saveOptions', saveOptions);
    setup;
    options.noiseType = strcat(['white gaussian 10^{-', num2str(x),'}']);
    actionReflectanceEstimationComparison;
    rmseCur = mean(rmse, 2);
end


function rmseCur = estimationIndependentNoise(x)
    dataset = 'saitama_v2';
    action = lower('ReflectanceEstimationSimple');
    skipLoading = false;
    showImages = false;
    saveImages = false;

    saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
    options = struct('tryReadData', false, 'dataset', dataset, 'action', action, ...
    'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
    'showImages', showImages, 'saveOptions', saveOptions);
    setup;
    options.noiseType = strcat(['independent 10^{-', num2str(x),'}']);
    actionReflectanceEstimationComparison;
    rmseCur = mean(rmse, 2);
end

function rmseCur = estimationSpatiallyAdaptive(x,y)
    dataset = 'saitama_v2';
    action = lower('ReflectanceEstimationSimple');
    skipLoading = false;
    showImages = false;
    saveImages = false;

    saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
    options = struct('tryReadData', false, 'dataset', dataset, 'action', action, ...
    'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
    'showImages', showImages, 'saveOptions', saveOptions);
    setup;
    options.noiseType = 'spatial';
    options.sigma1 = x;
    options.sigma2 = y;
    actionReflectanceEstimationComparison;
    rmseCur = mean(rmse, 2);
end