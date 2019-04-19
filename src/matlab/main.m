function [] = main(actions, showImages, saveImages, dataset, tryReadData)

%% Execute actions at bulk
%
% Usage: 
% main({'ReflectanceEstimationSimple', 'LBP'}, false, true, 'saitama_v2', false)
% 
% Input:
% 'action': string or cell array of strings, name of action to be performed 
%      actions = {
%          'SystemSpecs', 
%          'ReflectanceEstimationSystemComparison',
%          'ReflectanceEstimationMatrixComparison',
%          'ReflectanceEstimationMatrixSystemComparison',
%          'ReflectanceEstimationMatrixNoiseComparison',
%          'ReflectanceEstimationNoiseComparison',
%          'ReflectanceEstimationPreset', 
%          'ReflectanceEstimationSimple', 
%          'CreateSRGB',
%          'PCA',  
%          'LDA', 'LDA b'
%          'SVM', 
%          'KNN'
%          'LBP', 'LBP_RGB'
%          'EstimatedOverlap',
%          'ReflectanceEstimationOpposite', 
%          'PCALDA', 
%          'CountData'
%          'CompleteAnalysis',
%          'ClassificationPerformance'
%         };
% 'showImages' shows plots during execution (default: false)
% 'saveImages' save execution plots (default: true)
% 'dataset' string, name of the input dataset (default: 'saitama_v2')
%  datasets = {'color_sample', 'saitama', 'saitama_v2'}
% 'tryReadData' boolean, enables lada reading instead of loading (default:
% false) 
close all; %clc;

if (nargin < 2)
    showImages = false;
end
if (nargin < 3)
    saveImages = true;
end
if (nargin < 4)
    dataset = 'saitama_v7_min_region_e';
end
if (nargin < 5)
    tryReadData = false;
end

if ~iscell(actions) 
    actions = { actions };
else
    disp('Running batch execution');
end
    
for i = 1:numel(actions)
    %Set-up of options for running
    options =  setOpt([], dataset, actions{i}, showImages, saveImages, tryReadData);
    readData;

    %% main
    fprintf('Running action %s...\n', options.action);
    action = actions{i};
    if strcmp(action, 'SystemSpecs')
        actionCreateSystemPlots;
    elseif contains(action, 'ReflectanceEstimation')
        actionReflectanceEstimationComparison;
    elseif contains(action,'CreateSRGB')
        actionReconstructSRGB;
    elseif contains(action, 'PCA') || contains(action, 'LDA') || contains(action, 'DimRed')
        actionDimensionReduction;
    elseif contains(action, 'SVM') ||  contains(action, 'KNN')
       actionClassification;    
    elseif strcmp(action, 'ClassificationPerformance')
        actionClassificationPerformance;
    elseif contains(action, 'LBP')
        actionLBP;
    elseif contains(action, 'CompleteAnalysis')
        actionCompleteAnalysis;
    elseif contains(action, 'ApplyOnRGB')
        options.action = 'ReflectanceEstimationPreset_rgb';
        actionReflectanceEstimationComparison;
        options.action = 'lbp_rgb';
        actionLBP;
    else
        disp('Nothing to do here. Aborting...');
        return
    end
    
    tryReadData = false; %don't read again after the first time   
    logInfo(options);
    fprintf('Finished running action.\n');

end

disp('Finished execution')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [outputLog] = logInfo(options)
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
    for i = 1:width(optionsTab)
        outputLog = strjoin([outputLog, optionsTab.Properties.VariableNames(i), strrep(strrep(string(optionsTab(1, i).Variables), '..', ''), '\', ' '), '\n'], '     ');
    end
    outputLog = strjoin([outputLog, '-------------------------------------------\n\n\n\n\n\n\n']);
    logname = fullfile( '..', '..', 'logs', strcat( options.action, '_matlab.log'));
    fileID = fopen(fullfile('..', '..', 'logs', logname), 'a');
    fprintf(fileID, outputLog);
    fclose(fileID);

end