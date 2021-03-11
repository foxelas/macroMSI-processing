function [] = readWhite(dataDate, integrationTime, normByPixel, hasSmoothing, experiment, configuration, indirFolder)
%READWHITE reads necessary initialization info for white and black
%reference images
%
%   readWhite(dataDate, integrationTime, normByPixel, hasSmoothing,
%     experiment, configuration, indirFolder) prepares matfiles
%     with information about black and white reference images for
%     normalization of the HSI
%

close all;

if isempty(integrationTime)
    integrationTime = 1460;
end

%% Initialization
% normByPixel, hasSmoothing, experiment, configuration, indirFolder
%
initialization;

%% Filename
content = 'whiteReflectance';
hasFilter = ~(numel(strsplit(configuration, '_')) == 2); % eg configuration = 'singleLightClose_noFilter'
if ~hasFilter
    hasFilterSuffix = '_noFilter';
else 
    hasFilterSuffix = '';
end 
[filename, integrationTime] = getFilename(configuration, content, integrationTime);

%% Settings for normalization 1. Single Value, 2. PixelByPixel

if ~normByPixel
    setSetting('normFilename', mkNewDir(matdir, configuration, strcat(num2str(integrationTime), '_fullreflectance.mat')));
    setSetting('saveFolder', strcat(configuration, '_byAverage'));
else
    setSetting('normFilename', mkNewDir(matdir, configuration, strcat(num2str(integrationTime), '_fullreflectance_ByPixel.mat')));
    if hasSmoothing
        setSetting('saveFolder', strcat(configuration, '_', num2str(integrationTime), '_byPixel_withSmoothing'));
    else
        setSetting('saveFolder', strcat(configuration, '_', num2str(integrationTime), '_byPixel_noSmoothing'));
    end
end

normFilename = getSetting('normFilename');

%% Values for normalization
%For 99% reflectance
spectralData = loadH5Data(filename, configuration);
[~, ~, m] = size(spectralData);
wavelengths = getWavelengths(m);
whiteSize = size(spectralData);
dispImage = getDisplayImage(spectralData, 'rgb');
if ~normByPixel
    title('Select ROI for Normalization Spectrum');
    maskWhite = roipoly;
    fullReflectance = readHSI(spectralData, {maskWhite, maskWhite}, 'raw');
    save(normFilename, 'fullReflectance', 'maskWhite');
    fprintf('\nSaving white info at %s.\n\n', normFilename);
    setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('saveFolder'), 'whitePlot_'));
    plots(2, @plotSpectra, fullReflectance, wavelengths, '99%-white', 'Reflectance Spectrum of White Balance Sheet');
else
    xPoints = [100, 500, 900];
    if contains(configuration, 'singleLightFar')
        xPoints = xPoints + 400;
    end
    yPoints = [100, 500, 900];
    pointCount = numel(xPoints) * numel(yPoints);
    pointNames = cellfun(@(x) strcat('P', num2str(x)), num2cell(1:pointCount), 'UniformOutput', false);
    fullReflectanceByPixel = spectralData;
    pointSpectra = reshape(fullReflectanceByPixel(xPoints, yPoints, :), [pointCount, length(wavelengths)]);
    setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('saveFolder'), 'whitePlot_by_pixel_unsmoothened'));
    plots(2, @plotSpectra, pointSpectra, ...
        wavelengths, pointNames, 'Reflectance Spectrum of various pixels in White Balance Sheet');
    
    fullReflectanceByPixel = smoothImage(hasSmoothing, fullReflectanceByPixel);
    pointSpectra = reshape(fullReflectanceByPixel(xPoints, yPoints, :), [pointCount, length(wavelengths)]);
    setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('saveFolder'), 'whitePlot_by_pixel'));
    plots(2, @plotSpectra, pointSpectra, ...
        wavelengths, pointNames, 'Reflectance Spectrum of various pixels in White Balance Sheet');
    
    if sum(pointSpectra(2, :)' == squeeze(fullReflectanceByPixel(xPoints(2), yPoints(1), :))) ~= 401
        error('Not coinciding order of points');
    end
    
    plots(3, @plotPointsOnImage, dispImage, xPoints, yPoints, true);
    
    fprintf('\nSaving white info at %s.\n\n', normFilename);
    save(normFilename, 'fullReflectanceByPixel', '-v7.3');
end

%% Read Black
if nargin < 5
    indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'), 'h5');
else
    indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'), indirFolder, 'h5');
end
setSetting('datadir', indir);

filename = getFilename(strcat('noLight', hasFilterSuffix), 'capOn', integrationTime);

%For 99% reflectance
spectralData = loadH5Data(filename, configuration);
figure(3);
imagesc(spectralData(:, :, 100));
blackReflectance = spectralData;
blackReflectance = smoothImage(hasSmoothing, blackReflectance);

[x, y, z] = size(blackReflectance);
if (sum(size(blackReflectance) ~= whiteSize) > 0 && (size(blackReflectance, 1) < whiteSize(1)) && normByPixel)
    m = whiteSize(1);
    n = whiteSize(2);
    fullReflectanceByPixel = fullReflectanceByPixel(floor(m/2)-floor(x/2)+1:floor(m/2)+floor(x/2), floor(n/2)-floor(y/2)+1:floor(n/2)+floor(y/2), :);
end

savedir = mkNewDir(matdir, configuration, strcat(num2str(integrationTime), '_blackReflectance.mat'));
if ~normByPixel
    whiteMinusBlack = fullReflectance - reshape(blackReflectance, [x * y, z]);
    whiteMinusBlack = reshape(whiteMinusBlack, [x, y, z]);
    save(savedir, 'blackReflectance', 'whiteMinusBlack', '-v7.3');
else
    whiteMinusBlackByPixel = fullReflectanceByPixel - blackReflectance;
    save(savedir, 'blackReflectance', 'whiteMinusBlackByPixel', '-v7.3');
end
fprintf('\nSaved black info at %s.\n\n', savedir);

end

function result = smoothImage(hasSmoothing, target)
if hasSmoothing
    % smoothing the data with moving average filter
    windowSize = 10;
    b = (1 / windowSize) * ones(1, windowSize);
    a = 1;
    y = filter(b, a, target, [], 1);
    result = filter(b, a, y, [], 2);
else
    result = target;
end

end