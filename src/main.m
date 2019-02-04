
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
%          'ReflectanceEstimationPreset', 
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
    promt = strcat('Chooce an action:\n',...
        '1. Show System Specifications ', '2. Reflectance Estimation ', '3. Create sRGB', '4. Dimension Reduction ', ...
        '5. Classification ', '6. Extract LBP', '\n');
    selection = input(promt);

    if (selection == 1)
        actions = 'SystemPlots';
    elseif (selection == 2)
        actions = 'ReflectanceEstimationPreset';
    elseif (selection == 3)
        actions = 'CreatesRGB';
    elseif (selection == 4)
        actions = 'pca';
    elseif (selection == 5)
        actions = 'knn';
    elseif (selection == 6)
        actions = 'lbp';
    else
        error('Unsupported action.');
    end
       
end
if (nargin < 2)
    dataset = 'saitama_v2_min_region';
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

%% Set-up of options for running
saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
options = struct('tryReadData', tryReadData, 'dataset', dataset, 'action', actions{1}, ...
    'pixelValueSelectionMethod', 'extended', 'noiseType', 'sameForChannel', 'skipLoading', skipLoading, ...
    'showImages', showImages, 'saveOptions', saveOptions);
    
for i = 1:numel(actions)
    
    options = setOpt(options);
    readData;

    %% main
    tic;
    fprintf('Running action %s...\n', options.action);

    switch (options.action)
        case {'systemplots', 'systemspecs', 'SystemSpecs'}
            options.action = 'SystemSpecs';
            actionCreateSystemPlots;
            
        case lower('PrepareSmoothingMatrix') % needs to be performed before doing reflectance estimation 
            prepareSmoothingMatrix;

        case {'ReflectanceEstimationSystemComparison', ...
              'ReflectanceEstimationMatrixComparison', ...
              'ReflectanceEstimationMatrixSystemComparison', ...
              'ReflectanceEstimationNoiseComparison', ...
              'ReflectanceEstimationMatrixNoiseComparison', ...
              'ReflectanceEstimationPreset', ...
              'ReflectanceEstimationSimple'
              }
            actionReflectanceEstimationComparison;

        case {'CreatesRGB', 'reconstructsrgb'} % from the MSI reflectances   
            options.action = 'CreatesRGB';
            actionReconstructSRGB;

        case {'pca', 'lda', 'pcalda', 'dimred', 'pca b', 'lda b'}
            actionDimensionReduction;

        case {'knn'}
           options.action = 'KNNClassifier';
           actionClassification;
           
        case {'svm'}
           options.action = 'SVMClassifier';
           actionClassification;
           
        case {'lbp'}
            options.action = 'LBP';
            actionLBP;

        otherwise
            disp('Nothing to do here. Aborting...');
    end
    
    options.tryReadData = false; %don't read again after the first time
    
    fprintf('Finished running action.\n');

    elapsed = toc;

    logInfo(options, elapsed);
end

disp('Finished execution')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outputLog] = logInfo(options, runTime)
%LOGINFO produces a log of the current action run

    %% Make a struct without nested structs
    saveOptions = options.saveOptions;
    f = fieldnames(saveOptions);
    for i = 1:length(f)
        options.(f{i}) = saveOptions.(f{i});
    end
    options = rmfield(options, {'saveOptions'});

    %% Produce output log
    optionsTab = struct2table(options);
    outputLog = '-------------------------------------------\n';
    outputLog = strjoin([outputLog, 'Run on ', string(datetime), '\n'], ' ');
    outputLog = strjoin([outputLog, string(runTime), 's elapsed.', '\n'], ' ');
    for i = 1:width(optionsTab)
        outputLog = strjoin([outputLog, optionsTab.Properties.VariableNames(i), strrep(strrep(string(optionsTab(1, i).Variables), '..', ''), '\', ' '), '\n'], '     ');
    end
    outputLog = strjoin([outputLog, '-------------------------------------------\n\n\n\n\n\n\n']);
    logname = strcat(options.action, '_log.txt');
    fileID = fopen(fullfile('..', 'logs', logname), 'a');
    fprintf(fileID, outputLog);
    fclose(fileID);

end