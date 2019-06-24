function options = setOpt(options, dataset, showImages, saveImages, tryReadData)

if (isempty(options) && nargin == 5)
    saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
    options = struct('tryReadData', tryReadData, 'dataset', dataset, ...
        'pixelValueSelectionMethod', 'extended', 'noiseType', 'sameForChannel', ...
        'showImages', showImages, 'saveOptions', saveOptions);
end

if ~isfield(options, 'tryReadData')
    options.tryReadData = false;
end

if ~isfield(options, 'dataset')
    options.dataset = dataset;
end

if ~isfield(options, 'datadir')
    options.datadir = fullfile('..', '..', '..', '..', '..', '..', 'mspi', options.dataset);
end

if ~isfield(options, 'systemdir')
    options.systemdir = fullfile(getSource(), 'input', options.dataset);
end

if ~isfield(options, 'pixelValueSelectionMethod')
    options.pixelValueSelectionMethod = 'extended';
end

if ~isfield(options, 'saveOptions')
    options.saveOptions = struct('saveImages', true, 'saveInHQ', false, ...
        'savedir', fullfile(getSource(), 'output', options.dataset));
end

if ~isfield(options.saveOptions, 'savedir')
    options.saveOptions.savedir = fullfile(getSource(), 'output', options.dataset);
end

if ~isfield(options.saveOptions, 'saveInBW')
    options.saveOptions.saveInBW = false;
end

if ~isfield(options, 'showImages')
    options.showImages = false;
end

if ~isfield(options, 'luminanceCorrection')
    options.luminanceCorrection = 1;
end


validateattributes(options.luminanceCorrection,{'double'},{'positive','scalar', '<=', 1});

fprintf('Data directory is set to %s.\n', options.datadir);
fprintf('Save directory is set to %s.\n', options.saveOptions.savedir);

end