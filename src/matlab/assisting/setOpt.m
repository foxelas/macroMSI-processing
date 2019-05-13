function options = setOpt(options, dataset, action, showImages, saveImages, tryReadData)

if (isempty(options) && nargin == 6)
    saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
    options = struct('tryReadData', tryReadData, 'dataset', dataset, 'action', action, ...
        'pixelValueSelectionMethod', 'extended', 'noiseType', 'sameForChannel', ...
        'showImages', showImages, 'saveOptions', saveOptions);
end

if ~isfield(options, 'tryReadData')
    options.tryReadData = false;
end

if ~isfield(options, 'datadir')
    options.datadir = fullfile('..', '..', '..', '..', '..', '..', 'mspi', options.dataset);
end

if ~isfield(options, 'systemdir')
    options.systemdir = fullfile('..', '..', '..', 'input', options.dataset);
end

if ~isfield(options, 'pixelValueSelectionMethod')
    options.pixelValueSelectionMethod = 'extended';
end

if ~isfield(options, 'saveOptions')
    options.saveOptions = struct('saveImages', true, 'saveInHQ', false, ...
        'savedir', fullfile('..', '..', '..', 'output', options.dataset));
end

if ~isfield(options.saveOptions, 'savedir')
    options.saveOptions.savedir = fullfile('..', '..', '..','output', options.dataset);
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