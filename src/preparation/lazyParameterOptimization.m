close all; clc;

dataset = 'saitama_v3_min_region';
action = lower('ReflectanceEstimationSimple');
skipLoading = false;
showImages = true;
saveImages = false;

saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
options = struct('tryReadData', false, 'dataset', dataset, 'action', action, ...
'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
'showImages', showImages, 'saveOptions', saveOptions);

options = setOpt(options);
out = matfile(options.outName, 'Writable', true);
readData;

% %% For givenSNR
% snrRange = 5:30;
% errSNR = zeros(length(snrRange), 2);
% for i = 1:length(snrRange)
%     snr = snrRange(i);
%     options.noiseType = 'givenSNR';
%     options.smoothingMatrixMethod = 'Cor_All';
%     options.snr = snr;
%     actionReflectanceEstimationComparison;
%     errSNR(i,:) = [snr, minError];
% end
% minSnr = errSNR(errSNR(:,2) == min(errSNR(:,2)),1);

minSnr = 17;

%% For Adaptive
% gammaRange = 0.5:0.5:2;
% errGamma = zeros(length(gammaRange), 2);
% for i = 1:length(gammaRange)
%     gamma = gammaRange(i);
%     options.noiseType = 'givenSNR';
%     options.snr = minSnr; 
%     options.smoothingMatrixMethod = 'adaptive';
%     options.gamma = gamma;
%     actionReflectanceEstimationComparison;
%     errGamma(i,:) = [gamma, minError];
% end
minGamma = 2;


%% For spatially Adaptive 

sigma1Range = 0.001:0.05:0.2;
sigma2Range = 0.001:0.05:0.2;

errSigma = zeros(length(sigma1Range) * length(sigma2Range), 3);
for i = 1:length(sigma1Range)
    sigma1 = sigma1Range(i);
    for j = 1:length(sigma2Range)
        sigma2 = sigma2Range(j);
        options.noiseType = 'spatial';
        options.sigma1 = sigma1; 
        options.sigma2 = sigma2;
        options.smoothingMatrixMethod = 'Cor_All';
        actionReflectanceEstimationComparison;
        k = sub2ind([length(sigma1Range), length(sigma2Range)], i, j);
        errSigma(k,:) = [sigma1, sigma2, minError];

    end
end

