close all; clc;
dataset = 'saitama_v2_min_region';
action = lower('ReflectanceEstimationSimple');
skipLoading = false;
showImages = false;
saveImages = false;

saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
options = struct('tryReadData', false, 'dataset', dataset, 'action', action, ...
'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
'showImages', showImages, 'saveOptions', saveOptions);
setup;
fminsearchOptions = optimset('MaxIter', 5, 'PlotFcns',@optimplotfval);
fileID = fopen('..\logs\optimization.txt', 'a');

tic;
%% Optimize noise model "givenSNR"
[minSNR,minVal] = fminsearch(@estimationGivenSNR, 20, fminsearchOptions);
outputLog = sprintf('Given SNR:: MinRMSE=%.5f, SNR=%d\n', minVal, minSNR)
fprintf(fileID, outputLog);

%% Optimize noise model "adaptive"
[alpha,minVal] = fminsearch(@estimationAdaptive, 0.5, fminsearchOptions);
outputLog = sprintf('Adaptive Wiener:: MinRMSE=%.5f, Alpha=%.4f\n', minVal, alpha)
fprintf(fileID, outputLog);

% 
% %% Optimize noise model "White Gaussian"
% [minNoiseOrder,minVal] = fminsearch(@estimationWhiteGaussian, 2, fminsearchOptions);
% outputLog = sprintf('White gaussian:: MinRMSE=%.5f, NoiseOrder=10^(-%d)\n', minVal, minNoiseOrder)
% fprintf(fileID, outputLog);

%% Optimize noise model "SameForChannel noise"
[minSigma,minVal] = fminsearch(@estimationSameForChannel, 0.04, fminsearchOptions);
outputLog = sprintf('sigma n SNR:: MinRMSE=%.5f, sigma=%.5f\n', minVal, minSigma)
fprintf(fileID, outputLog);


%% Optimize noise model "DiffForChannel noise"
[sigma,minVal] = fminsearch(@estimationDiffForChannel, 0.03 * ones(7,1), fminsearchOptions);
outputLog = sprintf('Sigma noise:: MinRMSE=%.5f, Sigma=[%.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f]\n', minVal, sigma(1), sigma(2), sigma(3), sigma(4), sigma(5), sigma(6), sigma(7))
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

%% Optimize noise model "Spatial DiffForChannel"
[sigma,minVal] = fminsearch(@estimationSpatiallyAdaptive,  0.003 * ones(7,1), fminsearchOptions);
outputLog = sprintf('Spatially adaptive Wiener:: MinRMSE=%.5f, Sigma=[%.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f]\n', minVal, sigma(1), sigma(2), sigma(3), sigma(4), sigma(5), sigma(6), sigma(7))
fprintf(fileID, outputLog);

[sigma,minVal] = fminsearch(@estimationSpatiallyAdaptive, [0.05; 0.01], fminsearchOptions);
outputLog = sprintf('Spatially adaptive Wiener:: MinRMSE=%.5f, Sigma1=%.4f, Sigma2=%.4f\n', minVal, sigma(1), sigma(2))
fprintf(fileID, outputLog);


fclose(fileID);
t = toc;
fprintf('Action elapsed %.2f mins.\n', t / 60);
clear variables; 

function rmseCur = estimationGivenSNR(x)
    dataset = 'saitama_v2_min_region';
    action = lower('ReflectanceEstimationSimple');
    skipLoading = false;
    showImages = false;
    saveImages = false;

    saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
    options = struct('tryReadData', false, 'dataset', dataset, 'action', action, ...
    'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
    'showImages', showImages, 'saveOptions', saveOptions, 'smoothingMatrixMethod', 'Cor_Malignancy');
    readData;
    outName = setup(options);
    out = matfile(outName, 'Writable', true);
    options.noiseType = 'givenSNR';
    options.snr = x;
    actionReflectanceEstimationComparison;
    rmseCur = mean(rmse, 2);
end

function rmseCur = estimationAdaptive(x)
    dataset = 'saitama_v2_min_region';
    action = lower('ReflectanceEstimationSimple');
    skipLoading = false;
    showImages = false;
    saveImages = false;

    saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
    options = struct('tryReadData', false, 'dataset', dataset, 'action', action, ...
    'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
    'showImages', showImages, 'saveOptions', saveOptions);
    readData;
    outName = setup(options);
    out = matfile(outName, 'Writable', true);
    options.smoothingMatrixMethod = 'adaptive';
    options.noiseType = 'givenSNR';
    options.alpha = x;
    actionReflectanceEstimationComparison;
    rmseCur = mean(rmse, 2);
end

function rmseCur = estimationWhiteGaussian(x)
    dataset = 'saitama_v2_min_region';
    action = lower('ReflectanceEstimationSimple');
    skipLoading = false;
    showImages = false;
    saveImages = false;

    saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
    options = struct('tryReadData', false, 'dataset', dataset, 'action', action, ...
    'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
    'showImages', showImages, 'saveOptions', saveOptions, 'smoothingMatrixMethod', 'Cor_Malignancy');
    readData;
    outName = setup(options);
    out = matfile(outName, 'Writable', true);
    options.noiseType = strcat(['white gaussian 10^{-', num2str(x),'}']);
    actionReflectanceEstimationComparison;
    rmseCur = mean(rmse, 2);
end


function rmseCur = estimationDiffForChannel(x)
    dataset = 'saitama_v2_min_region';
    action = lower('ReflectanceEstimationSimple');
    skipLoading = false;
    showImages = false;
    saveImages = false;

    saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
    options = struct('tryReadData', false, 'dataset', dataset, 'action', action, ...
    'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
    'showImages', showImages, 'saveOptions', saveOptions);
    readData;
    outName = setup(options);
    out = matfile(outName, 'Writable', true);
    options.noiseType = 'diffForChannel';
    options.sigma = x';
    actionReflectanceEstimationComparison;
    rmseCur = mean(rmse, 2);
end

function rmseCur = estimationSameForChannel(x)
    dataset = 'saitama_v2_min_region';
    action = lower('ReflectanceEstimationSimple');
    skipLoading = false;
    showImages = false;
    saveImages = false;

    saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
    options = struct('tryReadData', false, 'dataset', dataset, 'action', action, ...
    'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
    'showImages', showImages, 'saveOptions', saveOptions, 'smoothingMatrixMethod', 'Cor_Malignancy');
    readData;
    outName = setup(options);
    out = matfile(outName, 'Writable', true);
    
    options.noiseType = 'sameForChannel';
    options.sigma = x;
    actionReflectanceEstimationComparison;
    rmseCur = mean(rmse, 2);
end

function rmseCur = estimationSpatiallyAdaptive(X)

    dataset = 'saitama_v2_min_region';
    action = lower('ReflectanceEstimationSimple');
    skipLoading = false;
    showImages = false;
    saveImages = false;

    saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
    options = struct('tryReadData', false, 'dataset', dataset, 'action', action, ...
    'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
    'showImages', showImages, 'saveOptions', saveOptions, 'smoothingMatrixMethod', 'Cor_Malignancy');
    setup;
    options.noiseType = 'spatial';
    if numel(X) == 2
        options.sigma1 = X(1);
        options.sigma2 = X(2);
    else
        options.variance = X';
    end
    
    actionReflectanceEstimationComparison;
    rmseCur = mean(rmse, 2);
end