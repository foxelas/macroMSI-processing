close all; clc;

dataset = 'saitama_v6_min_region_e';
action = 'ReflectanceEstimationSimple';
skipLoading = false;
showImages = false;
saveImages = false;

saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
options = struct('tryReadData', false, 'dataset', dataset, 'action', action, ...
'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
'showImages', showImages, 'saveOptions', saveOptions);

options = setOpt(options);
out = matfile(options.outName, 'Writable', true);
readData;

%% For givenSNR
 snrRange = 0:0.5:5;
 errTuning = struct('avgrmse', [], 'minrmse', [], 'maxrmse', [], ...
                'stdrmse',[], 'avgnmse', [], 'minnmse', [], ...
                'maxnmse', [], 'stdnmse', [], 'pixelValueSelectionMethod', [], ...
                'stdeSEM', [], ...
                'smoothingMatrixMethod', [], 'noiseType', []);
tuneValues = struct('SNR', []);
 for iter = 1:length(snrRange)
     snr = snrRange(iter);
     options.noiseType = 'givenSNR';
     options.smoothingMatrixMethod = 'Cor_Sample';
     options.noiseParam = snr;
     actionReflectanceEstimationComparison;
     errTuning(iter) = errors; 
     tuneValues(iter) = struct('SNR', snr);
 end

errTuning = catstruct(errTuning,tuneValues);
writetable(struct2table(errTuning),strcat('../../logs/spatially_adaptive_tuning_', options.dataset,'log.xlsx'), 'Sheet', 'SNR');
[~, minPos] = min([errTuning.avgrmse]);
minSnr = errTuning(minPos).SNR;

 %% For samesigma
sigmaRange = [10^(-8),10^(-7),  10^(-6),  10^(-5), 10^(-4),  10^(-3), 10^(-2), 10^(-1)];
errTuning = struct('avgrmse', [], 'minrmse', [], 'maxrmse', [], ...
                'stdrmse',[], 'avgnmse', [], 'minnmse', [], ...
                'maxnmse', [], 'stdnmse', [], 'pixelValueSelectionMethod', [], ...
                'stdeSEM', [], ...
                'smoothingMatrixMethod', [], 'noiseType', []);
tuneValues = struct('SameSigma', []);
 for iter = 1:length(sigmaRange)
     sigma = sigmaRange(iter);
     options.noiseType = 'sameForChannel';
     options.smoothingMatrixMethod = 'Cor_Sample';
     options.noiseParam = sigma;
     actionReflectanceEstimationComparison;
     errTuning(iter) = errors; 
     tuneValues(iter) = struct('SameSigma', sigma);
 end
errTuning = catstruct(errTuning,tuneValues);
writetable(struct2table(errTuning), strcat('../../logs/spatially_adaptive_tuning_', options.dataset,'log.xlsx'), 'Sheet', 'Same');
 
 
%% For Adaptive
 gammaRange = 0.5:0.5:2;
errTuning = struct('avgrmse', [], 'minrmse', [], 'maxrmse', [], ...
                'stdrmse',[], 'avgnmse', [], 'minnmse', [], ...
                'maxnmse', [], 'stdnmse', [], 'pixelValueSelectionMethod', [], ...
                'stdeSEM', [], ...
                'smoothingMatrixMethod', [], 'noiseType', []);
tuneValues = struct('Gamma', []);
 for iter = 1:length(gammaRange)
     gamma = gammaRange(iter);
     options.noiseType = 'givenSNR';
     options.noiseParam = minSnr;
     options.smoothingMatrixMethod = 'adaptive';
     options.gamma = gamma;
     actionReflectanceEstimationComparison;
     errTuning(iter) = errors; 
     tuneValues(iter) = struct('Gamma', gamma);
 end
errTuning = catstruct(errTuning,tuneValues);
writetable(struct2table(errTuning), strcat('../../logs/spatially_adaptive_tuning_', options.dataset,'log.xlsx'), 'Sheet', 'Adaptive');
 

%% For spatially Adaptive 

sigma2Range = [0.0001, 0.001, 0.01, 0.1];

errTuning = struct('avgrmse', [], 'minrmse', [], 'maxrmse', [], ...
                'stdrmse',[], 'avgnmse', [], 'minnmse', [], ...
                'maxnmse', [], 'stdnmse', [], 'pixelValueSelectionMethod', [], ...
                'stdeSEM', [], ...
                'smoothingMatrixMethod', [], 'noiseType', []);
sigmaValues = struct('sigma1', [], 'sigma2', []);
iter = 0;
for jj = 1:length(sigma2Range)
    sigma2 = sigma2Range(jj);
    sigma1Range = sigma2 * [0.01, 0.05, 0.1, 0.15, 0.5];
    for ii = 1:length(sigma1Range)
        sigma1= sigma1Range(ii);
       
        options.noiseType = 'spatial';
        options.noiseParam = [sigma1, sigma2]; 
        options.smoothingMatrixMethod = 'Cor_Sample';
        actionReflectanceEstimationComparison;
%         k = sub2ind([length(sigma1Range), length(sigma2Range)], i, j);
        iter = iter + 1;
        errTuning(iter) = errors; 
        sigmaValues(iter) = struct('sigma1', sigma1, 'sigma2', sigma2);
    end
end

options.noiseType = 'spatial Olympus';
actionReflectanceEstimationComparison;
iter = iter + 1;
errTuning(iter) = errors; 
sigmaValues(iter) = struct('sigma1', [], 'sigma2', []);

options.noiseType = 'fromOlympus';
actionReflectanceEstimationComparison;
iter = iter + 1;
errTuning(iter) = errors; 
sigmaValues(iter) = struct('sigma1', [], 'sigma2', []);

errTuning = catstruct(errTuning,sigmaValues);
writetable(struct2table(errTuning), strcat('../../logs/spatially_adaptive_tuning_', options.dataset,'log.xlsx'), 'Sheet', 'Spatial');


