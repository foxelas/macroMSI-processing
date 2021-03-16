function [] = readWhite(dataDate, integrationTime, experiment, configuration, indirFolder)
%READWHITE reads necessary initialization info for white and black
%reference images
%
%   readWhite(dataDate, integrationTime, configuration, indirFolder) prepares matfiles
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
setSetting('normFilename', mkNewDir(matdir, configuration, strcat(num2str(integrationTime), '_fullReflectance.mat')));
whiteFilename = getSetting('normFilename');

%% Values for normalization
%For 99% reflectance
spectralData = loadH5Data(filename, configuration);
[~, ~, m] = size(spectralData);
wavelengths = getWavelengths(m);
dispImage = getDisplayImage(spectralData, 'rgb');
setSetting('saveFolder', strcat(configuration, '_', num2str(integrationTime)));

%% uniSpectrum 
figure(1); 
imshow(dispImage);
title('Select ROI for Normalization Spectrum');
maskWhite = roipoly;
uniSpectrum = readHSI(spectralData, maskWhite, 'raw');
setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('saveFolder'), 'whitePlot_unispectrum'));
plots(2, @plotSpectra, uniSpectrum, wavelengths, '99%-white', 'Reflectance Spectrum of White Balance Sheet');
    
%% byPixel 
xPoints = [100, 500, 900];
if contains(configuration, 'singleLightFar')
    xPoints = xPoints + 400;
end
yPoints = [100, 500, 900];
pointCount = numel(xPoints) * numel(yPoints);
pointNames = cellfun(@(x) strcat('P', num2str(x)), num2cell(1:pointCount), 'UniformOutput', false);
fullReflectanceByPixel = spectralData;
pointSpectra = reshape(fullReflectanceByPixel(xPoints, yPoints, :), [pointCount, length(wavelengths)]);
setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('saveFolder'), 'whitePlot_byPixel'));
plots(3, @plotSpectra, pointSpectra, ...
    wavelengths, pointNames, 'Reflectance Spectrum of various pixels in White Balance Sheet');

if sum(pointSpectra(2, :)' == squeeze(fullReflectanceByPixel(xPoints(2), yPoints(1), :))) ~= 401
    error('Not coinciding order of points');
end

plots(4, @plotPointsOnImage, dispImage, xPoints, yPoints, true);

fprintf('\nSaving white info at %s.\n\n', whiteFilename);
save(whiteFilename, 'fullReflectanceByPixel', 'uniSpectrum', 'maskWhite', '-v7.3');


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
dispImage = getDisplayImage(spectralData, 'rgb');
figure(5); imshow(dispImage);
setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('saveFolder'), 'blackImage'));
savePlot(5);

blackReflectance = spectralData;
savedir = mkNewDir(matdir, configuration, strcat(num2str(integrationTime), '_blackReflectance.mat'));
save(savedir, 'blackReflectance', '-v7.3');
fprintf('\nSaved black info at %s.\n\n', savedir);

end
