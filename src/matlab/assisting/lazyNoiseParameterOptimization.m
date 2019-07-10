function [] = lazyNoiseParameterOptimization(dataset)
%%lazyNoiseParameterOptimization perform lazy parameter optimization for 
% the noise model in reflectance estimation 

    action = 'ReflectanceEstimationSimple';
    showImages = false;
    saveImages = false;
    saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
    
    options = struct('tryReadData', false, 'dataset', dataset, 'action', action, ...
    'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', ...
    'showImages', showImages, 'saveOptions', saveOptions);
    
    minValue = optimizeNoiseParameters('snr', options);
    optimizeNoiseParameters('sameForChannel', options);
    optimizeNoiseParameters('adaptive', options, minValue);
    optimizeNoiseParameters('spatial', options);

end

function [minValue] = optimizeNoiseParameters(noiseModel, options, minSnr)

minValue = [];
if nargin < 3 
    minSnr = 0;
end

options = setOpt(options);
readData;
    
[tuneValues, tuneRange] = getTuneValueStruct(noiseModel);

 errTuning = struct('avgrmse', [], 'minrmse', [], 'maxrmse', [], ...
            'stdrmse',[], 'avgnmse', [], 'minnmse', [], ...
            'maxnmse', [], 'stdnmse', [], 'pixelValueSelectionMethod', [], ...
            'stdeSEM', [], 'smoothingMatrixMethod', [], 'noiseType', []);
            
if contains(noiseModel, 'snr')
    for iter = 1:length(tuneRange)
         snr = tuneRange(iter);
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
    minValue = errTuning(minPos).SNR;    
    
elseif contains(noiseModel, 'sameForChannel')
     for iter = 1:length(tuneRange)
         sigma = tuneRange(iter);
         options.noiseType = 'sameForChannel';
         options.smoothingMatrixMethod = 'Cor_Sample';
         options.noiseParam = sigma;
         actionReflectanceEstimationComparison;
         errTuning(iter) = errors; 
         tuneValues(iter) = struct('SameSigma', sigma);
     end
    errTuning = catstruct(errTuning,tuneValues);
    writetable(struct2table(errTuning), strcat('../../logs/spatially_adaptive_tuning_', options.dataset,'log.xlsx'), 'Sheet', 'Same');


elseif contains(noiseModel, 'adaptive')
     for iter = 1:length(tuneRange)
         gamma = tuneRange(iter);
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

elseif contains(noiseModel, 'spatial')
    iter = 0;
    for jj = 1:length(tuneRange)
        sigma2 = tuneRange(jj);
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
            tuneValues(iter) = struct('sigma1', sigma1, 'sigma2', sigma2);
        end
    end

    options.noiseType = 'spatial Olympus';
    actionReflectanceEstimationComparison;
    iter = iter + 1;
    errTuning(iter) = errors; 
    tuneValues(iter) = struct('sigma1', [], 'sigma2', []);

    options.noiseType = 'fromOlympus';
    actionReflectanceEstimationComparison;
    iter = iter + 1;
    errTuning(iter) = errors; 
    tuneValues(iter) = struct('sigma1', [], 'sigma2', []);

    errTuning = catstruct(errTuning,tuneValues);
    writetable(struct2table(errTuning), strcat('../../logs/spatially_adaptive_tuning_', options.dataset,'log.xlsx'), 'Sheet', 'Spatial');

else 
    error('Unknown noiseModel')
end
               
end

function [tuneValues, tuneRange] = getTuneValueStruct(noiseModel)

if contains(noiseModel, 'snr')
    tuneValues = struct('SNR', []);
    tuneRange = 0:0.5:5;
elseif contains(noiseModel, 'sameForChannel')
    tuneValues = struct('SameSigma', []);
    tuneRange = [10^(-8),10^(-7),  10^(-6),  10^(-5), 10^(-4),  10^(-3), ...
        10^(-2), 10^(-1)];
elseif contains(noiseModel, 'adaptive')
    tuneValues = struct('Gamma', []);
    tuneRange = 0.5:0.5:2;
elseif contains(noiseModel, 'spatial')
    tuneValues = struct('sigma1', [], 'sigma2', []);
    tuneRange = [0.0001, 0.001, 0.01, 0.1];
else 
    error('Unknown noiseModel')
end
               
end