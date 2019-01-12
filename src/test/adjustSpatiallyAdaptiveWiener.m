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
tic;

fileID = fopen('..\logs\log.txt', 'a');
minVal = 2;
options.noiseType = 'givenSNR';
for snr = 35:0.1:37
    options.snr = snr;
    actionReflectanceEstimationComparison;
    minCur = mean(rmse, 2);
    if minCur < minVal
        minVal = minCur;
        minSNR = snr;
    end
end
outputLog = sprintf('Given SNR:: MinRMSE=%.5f, SNR=%d\n', minVal, minSNR)
fprintf(fileID, outputLog);
clear('rmse');

minVal = 2;
for noiseOrder = 2:6
    options.noiseType = strcat(['white gaussian 10^{-', num2str(noiseOrder),'}']);
    actionReflectanceEstimationComparison;
    minCur = mean(rmse, 2);
    if minCur < minVal
        minVal = minCur;
        minNoiseOrder = noiseOrder;
    end
end
outputLog = sprintf('White gaussian:: MinRMSE=%.5f, NoiseOrder=10^(-%d)\n', minVal, minNoiseOrder)
fprintf(fileID, outputLog);
clear('rmse');

minVal = 2;
for noiseOrder = 3:8
    options.noiseType = strcat(['independent 10^{-', num2str(noiseOrder),'}']);
    actionReflectanceEstimationComparison;
    minCur = mean(rmse, 2);
    if minCur < minVal
        minVal = minCur;
        minNoiseOrder = noiseOrder;
    end
end
outputLog = sprintf('Independent noise:: MinRMSE=%.5f, NoiseOrder=10^(-%d)\n', minVal, minNoiseOrder)
fprintf(fileID, outputLog);
clear('rmse');

options.noiseType = 'fromOlympus';
actionReflectanceEstimationComparison;
outputLog = sprintf('Simple Wiener:: MinRMSE=%.5f, Noise from Olympus\n', mean(rmse, 2))
fprintf(fileID, outputLog);
clear('rmse');

options.noiseType = 'spatial';
actionReflectanceEstimationComparison;
outputLog = sprintf('Spatially adaptive Wiener:: MinRMSE=%.5f, Noise from Olympus\n', mean(rmse, 2))
fprintf(fileID, outputLog);
clear('rmse');

minVal = 2;
options.noiseType = 'spatial';
for sigma1 = 10^(-4):(5*10^(-4)):0.03
    for sigma2 = 10^(-4):(5*10^(-4)):0.3
        options.sigma1 = sigma1;
        options.sigma2 = sigma2;
        actionReflectanceEstimationComparison;
        minCur = mean(rmse, 2);
        if minCur < minVal 
            minVal = minCur;
            sigma1Min = sigma1;
            sigma2Min = sigma2;
        end
    end
end
outputLog = sprintf('Spatially adaptive Wiener:: MinRMSE=%.5f, Sigma1=%.4f, Sigma2=%.4f\n', minVal, sigma1Min, sigma2Min)
fprintf(fileID, outputLog);

fclose(fileID);

t = toc;
fprintf('Action elapsed %.2f mins.\n', t / 60);
clear variables; 

