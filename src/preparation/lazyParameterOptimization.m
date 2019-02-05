close all; clc;

dataset = 'saitama_v5_min_square';
action = lower('ReflectanceEstimationSimple');
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

% %% For givenSNR
%  snrRange = 5:5:30;
%  errSNR = zeros(length(snrRange), 2);
%  for i = 1:length(snrRange)
%      snr = snrRange(i);
%      options.noiseType = 'givenSNR';
%      options.smoothingMatrixMethod = 'Cor_Sample';
%      options.noiseParam = snr;
%      actionReflectanceEstimationComparison;
%      errSNR(i,:) = [snr, minError];
%  end
%  minSnr = errSNR(errSNR(:,2) == min(errSNR(:,2)),1)
% 
%  %% For samesigma
%  sigmaRange = [10^(-8),10^(-7),  10^(-6),  10^(-5), 10^(-4),  10^(-3), 10^(-2), 10^(-1)];
%  errSigma = zeros(length(sigmaRange), 2);
%  for i = 1:length(sigmaRange)
%      sigma = sigmaRange(i);
%      options.noiseType = 'sameForChannel';
%      options.smoothingMatrixMethod = 'Cor_Sample';
%      options.noiseParam = sigma;
%      actionReflectanceEstimationComparison;
%      errSigma(i,:) = [sigma, minError];
%  end
%  minSigma = errSigma(errSigma(:,2) == min(errSigma(:,2)),1)
 
 
%% For Adaptive
 gammaRange = 0.5:0.5:2;
 errGamma = zeros(length(gammaRange), 2);
 for i = 1:length(gammaRange)
     gamma = gammaRange(i);
     options.noiseType = 'givenSNR';
     options.noiseParam = minSnr;
     options.smoothingMatrixMethod = 'adaptive';
     options.gamma = gamma;
     actionReflectanceEstimationComparison;
     errGamma(i,:) = [gamma, minError];
 end
minGamma = errGamma(errGamma(:,2) == min(errGamma(:,2)),1)



%% For spatially Adaptive 

sigma1Range = 0.001:0.005:0.2;
sigma2Range = 0.01:0.01:0.2;

errSigma = zeros(length(sigma1Range) * length(sigma2Range), 3);
for i = 1:length(sigma1Range)
    sigma1 = sigma1Range(i);
    for j = 1:length(sigma2Range)
        sigma2 = sigma2Range(j);
        options.noiseType = 'spatial';
        options.noiseParam = [sigma1, sigma2]; 
        options.smoothingMatrixMethod = 'Cor_Sample';
        actionReflectanceEstimationComparison;
        k = sub2ind([length(sigma1Range), length(sigma2Range)], i, j);
        errSigma(k,:) = [sigma1, sigma2, minError];
    end
end
[~, idMin] = errSigma(:,3);
minSigma12 = errSigma(idMin,1:2);


