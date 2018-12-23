function outputLog  = batchExecution( dataset, skipLoading, showImages, saveImages )
%% Execute actions at bulk
% 
% Input 
% 
% dataset: the input dataset (default: 'saitama_v2') 
% skipLoading: see below (default: false)
% showImages: see below (default: false)
% saveImages: see below (default: true)
% 
% Options structure:
%
% 'tryReadData' boolean, enables lada reading instead of loading (default: false)
% 'dataset' string, name of the input dataset (default: 'saitama_v2')
%     datasets = {'color_sample', 'saitama', 'saitama_v2'}
% 'action' string, name of action to be performed (default: 'ConfPlots')
%      actions = {'ConfPlots', 'FixMinimumErrorPoints', 'ReflectanceEstimationSystemComparison', 'ReflectanceEstimationMatrixComparison', ...
%           'ReflectanceEstimationMatrixSystemComparison', 'ReflectanceEstimationNoiseComparison', 'ReflectanceEstimationSimple', 'CreateSRGB', ...
%           'PCA' , 'LDA', 'Nearest', 'MeasuredOverlap', 'EstimatedOverlap', 'ReflectanceEstimationOpposite', 'PCA-LDA', 'CountData'};
% 'classification' string, name of the classification method to be applied (default: 'knn')
%     classification = { 'knn' };
% 'smoothingMatrixMethod' string, name of the smoothing matrix of wiener estimation (default: 'corr same fixing all spectra')
%      smoothingMatrixMethods = {'markovian', 'corr all spectra', 'corr macbeth spectra', 'corr sample spectra', ...
%         'corr same type all spectra', 'corr same type sample spectra', 'corr same fixing all spectra', 'corr same fixing sample spectra'};
% 'pixelValueSelectionMethods' string, name of the pixel value selection method to reduce the 4D raw input to 3D MSI (default: 'extended')
%      pixelValueSelectionMethods = {'green', 'rms', 'adjusted', 'extended'}; 
% 'noiseModel' string, name of the node model to be used (default: 'givenSNR')
%      noiseModels = {'independent', 'fromOlympus', 'dependent', 'givenSNR', 'none'};
% 'snr' double, noise Signal-To-Noise Ratio (default: 25) 
% 'skipLoading' boolean, disables re-loading of workspace input data to save time (default: true)
% 'showImages' boolean, enables figure preview while running (default: false)
% 'saveOptions' struct, manages output saving options
%     saveOptions = struct('saveImages', true, 'saveInHQ', false)
close all; clc;

if (nargin < 1) 
    dataset = 'saitama_v2';
end
if (nargin < 2)
    skipLoading = false;
end
if (nargin < 3)
    showImages = false;
end
if (nargin < 4)
    saveImages = true;
end

saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
options = struct('tryReadData', false, 'dataset', dataset, 'action', 'ReflectanceEstimationMatrixSystemComparison', ...
    'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
    'showImages', showImages, 'saveOptions', saveOptions);

disp('Running batch execution')
main 
% outputLog is produced inside main;
disp('Finished batch execution')

end

% options.action = 'ReflectanceEstimationSimple';
% main
% options.action = 'ReflectanceEstimationSystemComparison';
% main
% options.action = 'ReflectanceEstimationNoiseComparison';
% main
% options.action = 'ReflectanceEstimationMatrixComparison';
% main
% options.action = 'ReflectanceEstimationMatrixSystemComparison';
% main
% options.action = 'LDA';
% main
% options.action = 'PCA';
% main
% options.action = 'PCA-LDA';
% main
% options.action = 'Nearest';
% main
% options.action = 'ReflectanceEstimationOpposite';
% main
% options.action = 'CreateSRGB';
% main

