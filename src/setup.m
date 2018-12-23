%% SET action OPTIONS
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
%
% options = struct(...)
% e.g. options = struct('tryReadData', false, 'dataset', dataset, 'action', 'ReflectanceEstimationSimple', ...
%     'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'skipLoading', skipLoading, ...
%     'showImages', showImages, 'saveOptions', saveOptions);

if ~isfield(options, 'tryReadData')
    options.tryReadData = false;
end 

if ~isfield(options, 'dataset')
    error('No data selected');
end 

if ~isfield(options, 'datadir')
    options.datadir = fullfile('..', '..', '..', '..', '..', 'mspi', dataset);
end
options.systemdir = fullfile('..', '..', 'input', options.dataset);

if ~isfield(options, 'action')
    options.action = 'ConfPlots';
end

if ~isfield(options, 'smoothingMatrixMethod')
    options.smoothingMatrixMethod = 'markovian';
end

if ~isfield(options, 'noiseType')
    options.noiseType = 'givenSNR';
end

if ~isfield(options, 'pixelValueSelectionMethod')
    options.pixelValueSelectionMethod = 'extended';
end

if ~isfield(options, 'saveOptions')
    options.saveOptions = struct('saveImages', true, 'saveInHQ', false, ...
        'savedir', fullfile('..', '..', 'output', options.dataset));
end

if ~isfield(options.saveOptions, 'savedir')
    options.saveOptions.savedir = fullfile('..', '..', 'output', options.dataset);  
end

if ~isfield(options, 'skipLoading')
    options.skipLoading = true;
end

if ~isfield(options, 'showImages')
    options.showImages = false;
end

fprintf('Data directory is set to %s.\n', options.datadir);
fprintf('Save directory is set to %s.\n', options.saveOptions.savedir);

pixelValueSelectionMethods = {'green', 'rms', 'adjusted', 'extended', 'rgb'};

smoothingMatrixMethods = {'markovian', 'corr all spectra', 'corr macbeth spectra', 'corr sample spectra', ...
    'corr same type all spectra', 'corr same type sample spectra', 'corr same fixing all spectra', 'corr same fixing sample spectra'};
nms = {'none', 'white gaussian 10^{-3}', 'white gaussian 10^{-5}', 'independent 10^{-3}', 'independent 10^{-5}', 'givenSNR 10dB', 'givenSNR 15dB', 'fromOlympus'}; 

fn = fullfile(options.saveOptions.savedir, options.action);
if ~exist(fn, 'dir')
    mkdir(fn);
    addpath(fn);
end

outName = generateName(options, 'matfileout');
out = matfile(outName,'Writable',true);

% SET action OPTIONS ends

%% LOAD DATA & INITIALIZE
if ~options.skipLoading
    disp('Initialization.')  
    load(fullfile(options.systemdir, 'system.mat')); % camera system parameters 
    load(fullfile(options.systemdir, 'data.mat')); % image data
    load(fullfile(options.systemdir,'ID.mat'));    % image data id and info struct
    
    msiN = length(ID);
    if (options.tryReadData) 
        fprintf('Reading dataset...\n\n');
        %readAndSaveSpectrum(options, ID, data);
        readAndSaveMSI(options, ID, data, 3, 3, bands); 
    end
%     load(strcat(options.systemdir,'coeff.mat'));
    in = matfile(generateName(options, 'matfilein')');
    MSIStruct = in.MSIStruct;
    whiteStruct = in.WhiteMSIStruct;
    darkStruct = in.DarkMSIStruct;
    measuredSpectrumStruct = in.MeasuredSpectrumStruct;
    name = options.dataset;
    fprintf('Finished initialization.\n\n')  
    fprintf('Let''s get down to business :muscle:\n\n');
end
% LOAD DATA & INITIALIZE ends