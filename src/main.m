function [] = main(actions, dataset, skipLoading, showImages, saveImages, tryReadData)

%% Execute actions at bulk
%
% Usage: 
% main('ReflectanceEstimationSimple', 'saitama_v2', false, true, false, false)
% 
% Input
%
% actions: string or cell array of strings
% dataset: the input dataset (default: 'saitama_v2') skipLoading: see below
% (default: false) showImages: see below (default: false) saveImages: see
% below (default: true)
%
% Options structure:
%
% 'tryReadData' boolean, enables lada reading instead of loading (default:
% false) 'dataset' string, name of the input dataset (default:
% 'saitama_v2')
%     datasets = {'color_sample', 'saitama', 'saitama_v2'}
% 'action' string, name of action to be performed (default: 'ConfPlots')
%      actions = {
%          'ConfPlots', 
%          'FixMinimumErrorPoints',
%          'ReflectanceEstimationSystemComparison',
%          'ReflectanceEstimationMatrixComparison',
%          'ReflectanceEstimationMatrixSystemComparison',
%          'ReflectanceEstimationMatrixNoiseComparison',
%          'ReflectanceEstimationNoiseComparison',
%          'ReflectanceEstimationSimple', 
%          'CreateSRGB',
%          'PCA' , 
%          'LDA',
%          'classification', 
%          'segmentation'
%          'MeasuredOverlap', 
%          'EstimatedOverlap',
%          'ReflectanceEstimationOpposite', 
%          'PCALDA', 
%          'CountData'
%         };
% 'classification' string, name of the classification method to be applied
% (default: 'knn')
%     classification = { 'knn' };
% 'smoothingMatrixMethod' string, name of the smoothing matrix of wiener
% estimation (default: 'corr same fixing all spectra')
%      smoothingMatrixMethods = {'markovian', 'Cor_All', 'Cor_Macbeth', 'Cor_Malignancy',...
%                                 'Cor_Fixing', 'Cor_MalignancyFixing', 
%                                 'Cor_Sample', 'Cor_SampleMalignancyFixing'};
% 'pixelValueSelectionMethods' string, name of the pixel value selection
% method to reduce the 4D raw input to 3D MSI (default: 'extended')
%      pixelValueSelectionMethods = {'green', 'rms', 'adjusted',
%      'extended'};
% 'noiseModel' string, name of the node model to be used (default:
% 'givenSNR')
%      noiseModels = {'independent', 'fromOlympus', 'dependent',
%      'givenSNR', 'none'};
% 'snr' double, noise Signal-To-Noise Ratio (default: 25) 'skipLoading'
% boolean, disables re-loading of workspace input data to save time
% (default: true) 'showImages' boolean, enables figure preview while
% running (default: false) 'saveOptions' struct, manages output saving
% options
%     saveOptions = struct('saveImages', true, 'saveInHQ', false)

close all; clc;
if (nargin < 1)
    actions = 'segmentation';
end
if (nargin < 2)
    dataset = 'saitama_v2';
end
if (nargin < 3)
    skipLoading = false;
end
if (nargin < 4)
    showImages = false;
end
if (nargin < 5)
    saveImages = true;
end
if (nargin < 6)
    tryReadData = false;
end

if ~iscell(actions) 
    actions = { actions };
else
    disp('Running batch execution');
end

for i = 1:numel(actions)
    
    action = lower(actions{i});
    
    %% Set-up of options for running
    saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
    options = struct('tryReadData', tryReadData, 'dataset', dataset, 'action', action, ...
        'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
        'showImages', showImages, 'saveOptions', saveOptions);
    setup;

    %% main
    tic;
    fprintf('Running action %s...\n', options.action);

    switch (options.action)
        case {'confplots', 'configureplots'}
            actionCreateConfPlots;

        case 'countdata'
            actionCountData;
            
        case lower('prepareSmoothingMatrix')  % needs to be performed before doing reflectance estimation 
            prepareSmoothingMatrix;

        case lower('FixMinimumErrorPoints')
            actionFindMinimumRefErrorPoints;

        case {lower('ReflectanceEstimationSystemComparison'), ...
              lower('ReflectanceEstimationMatrixComparison'), ...
              lower('ReflectanceEstimationMatrixSystemComparison'), ...
              lower('ReflectanceEstimationNoiseComparison'), ...
              lower('ReflectanceEstimationMatrixNoiseComparison'), ...
              lower('ReflectanceEstimationSimple')}
            actionReflectanceEstimationComparison;

        case {'createsrgb', 'reconstructsrgb'} % from the MSI reflectances       
            actionReconstructSRGB;

        case {'pca', 'lda', 'pcalda', 'dimred'}
            actionDimensionReduction;

        case {'classification', 'class', 'nearest', 'knn'}
           actionClassification;

        case {'measuredoverlap', 'estimatedoverlap'}
           actionOverlap;

        otherwise
            disp('Nothing to do here. Aborting...');
    end

    fprintf('Finished running action %s.\n', options.action);

    elapsed = toc;

    logInfo(options, elapsed);
end

disp('Finished execution')

end

