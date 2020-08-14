function [melMap, hbMap] = getMap(msi, mapType, mask, id, msiType)
%     GETMAP returns an optical map showing chromophore components of the macropathology image
%
%     Usage:
%     [melMap, hbMap] = getMap(msi)
%     [melMap, hbMap] = getMap(msi, 'ding')
%
%     Available mapType values
%     ding: uses optical density
%     vasefi: uses absorption slope
%     ours: uses estimated reflectance and absorption slope

if nargin < 2
    mapType = 'ding';
end

[channels, height, width] = size(msi);

if nargin < 3
    mask = ones(height, width);
end

if nargin < 5
    msiType = 'adjusted';
end

% reference = getReference(getSetting('systemdir'), height, width);
reference = getReferenceFromMacbeth(height, width);
reference = raw2msi(reference, msiType);

savedir = getSetting('savedir');
mapdir = getSetting('map');
% setSetting('plotName', fullfile(savedir, mapdir, 'msi.png'));
% plotFunWrapper(2, @plotMSI, msi);
mask3d = permute(repmat(mask, [1, 1, size(msi, 1)]), [3, 1, 2]);
msi(~mask3d) = nan;

normalizedReflectance = msi ./ reference;
% setSetting('plotName', fullfile(savedir, mapdir, 'normalizedMsi.png'));
% plotFunWrapper(3, @plotMSI, normalizedReflectance);

logIm = log10(normalizedReflectance);
opticalDensity = logIm;
absorption = -logIm;
% setSetting('plotName', fullfile(savedir, mapdir, 'absoprtion.png'));
% plotFunWrapper(4, @plotMSI, absorption);
switch mapType
    case 'ding'
        melMap = squeeze(opticalDensity(7, :, :));
        melBarTitle = 'Optical Density of Melanin';
        melSaveName = 'DingMel';
        melMapLimits = [-1.5, 0.5];
        hbMap = squeeze(opticalDensity(5, :, :)) - 1.15 .* squeeze(opticalDensity(7, :, :));
        hbBarTitle = 'Optical Density of Hemoglobin';
        hbMapSaveName = 'DingHb';
        hbMapLimits = [-5, 1];

    case 'vasefi'
        fc = [450, 465, 505, 525, 575, 605, 630]';
        range = [605, 630];
        melMap = estimateSlope(absorption, range, fc);
        setSetting('plotName', fullfile(savedir, mapdir, 'scaledMap.png'));
        melBarTitle = 'Melanin Absorbance Slope (a.u.)';
        melSaveName = 'VasefiMel';
        melMapLimits = [-0.06, 0.01];
        range = [505, 575];
        hbMap = estimateSlope(absorption, range, fc);
        setSetting('plotName', fullfile(savedir, mapdir, 'scaledMap.png'));
        hbBarTitle = 'Hemoglobin Absorbance Slope (a.u.)';
        hbMapSaveName = 'VasefiHb';
        hbMapLimits = [-0.06, 0.06];

    case 'diebele'
        melMap = squeeze(opticalDensity(7, :, :));
        melBarTitle = 'Melanin Index';
        melSaveName = 'DiebeleMel';
        melMapLimits = [-1.5, 0.5];
        hbMap = squeeze(opticalDensity(7, :, :)) - squeeze(opticalDensity(5, :, :));
        hbBarTitle = 'Erythema Index';
        hbMapSaveName = 'DiebeleHb';
        hbMapLimits = [-1, 5];

    case 'kapsokalyvas'
        %totalRed = squeeze(sum(squeeze(raw(:,:,:,3)), 1));
        %totalGreen = squeeze(sum(squeeze(raw(:,:,:,2)), 1));
        totalRed = squeeze(msi(7, :, :));
        totalGreen = squeeze(msi(5, :, :));
        totalBlue = squeeze(msi(1, :, :));
        melMap = (totalRed - min(totalRed(:))) ./ (max(totalRed(:)) - min(totalRed(:)));
        melBarTitle = 'Melanin Contrast';
        melMap = (totalBlue - totalGreen) ./ (totalBlue + totalGreen);
        melBarTitle = 'Superficial Melanin';
        melSaveName = 'KapsokalyvasMel';
        melMapLimits = [-0.6, 1];
        hbMap = (totalGreen - totalRed) ./ (totalGreen + totalRed);
        hbBarTitle = 'Hemoglobin Contrast';
        hbMapSaveName = 'KapsokalyvasHb';
        hbMapLimits = [-1, 0.8];

    case 'kuzmina'
        load(getSetting('channelCorrection'), 'redCoeff', 'greenCoeff', 'blueCoeff');
        calibCoeff = 1;
        totalRed = squeeze(msi(7, :, :));
        totalGreen = squeeze(msi(5, :, :));
        totalBlue = squeeze(msi(1, :, :));
        calibCoeff = blueCoeff / redCoeff;
        melMap = calibCoeff * totalBlue ./ totalRed;
        melBarTitle = 'Melanin Index';
        melSaveName = 'KuzminaMel';
        melMapLimits = [0, 25];
        calibCoeff = redCoeff / greenCoeff;
        hbMap = calibCoeff * totalRed ./ totalGreen;
        hbBarTitle = 'Hemoglobin Index';
        hbMapSaveName = 'KuzminaHb';
        hbMapLimits = [];

    case 'ours'
        [absorption, wavelength] = getAbsorptionByEstimatedReflectance(msi, mask, id);
        range = [605, 630];
        load(fullfile(getSetting('systemdir'), 'system.mat'), 'wavelength');
        melMap = estimateSlope(absorption, range, wavelength);
        melBarTitle = 'Melanin Absorbance Slope (a.u.)';
        melSaveName = 'oursMel';
        melMapLimits = [-1.49, -1.472] .* 0.001;
        range = [505, 575];
        hbMap = estimateSlope(absorption, range, wavelength);
        hbBarTitle = 'Hemobglobin Absorbance Slope (a.u.)';
        hbMapSaveName = 'oursHb';
        hbMapLimits = [1.41, 1.50] .* 0.001;

    otherwise
        disp('Unsupported type')
end
melMap = removeInf(melMap);
hbMap = removeInf(hbMap);
if isstruct(id)
    id = id.Index;
end
setSetting('plotName', fullfile(savedir, mapdir, strcat(melSaveName, '_', num2str(id), '_', 'ScaledMap.png')));
plotFunWrapper(1, @plotMap, melMap, mask, [], false, melBarTitle, melMapLimits);
setSetting('plotName', fullfile(savedir, mapdir, strcat(hbMapSaveName, '_', num2str(id), '_', 'ScaledMap.png')));
plotFunWrapper(2, @plotMap, hbMap, mask, [], false, hbBarTitle, hbMapLimits);

end

function slopes = estimateSlope(image, range, fc)

rangeStart = range(1);
rangeEnd = range(2);
[b, m, n] = size(image);
condition = fc >= rangeStart & fc <= rangeEnd;
x = fc(condition);
columnImage = reshape(image, b, m*n);
y = columnImage(condition, :);
X = [ones(length(x), 1), x];
b = X \ y;
slopes = b(2, :);
slopes = reshape(slopes, m, n);

end

function [absorption, wavelength] = getAbsorptionByEstimatedReflectance(I, mask, id)
setSetting('pixelValueSelectionMethod', 'adjusted');
setSetting('smoothingMatrixMethod', 'Cor_All');
setSetting('noiseType', 'fromOlympus');
estimatedReflectance = estimateReflectance(I, mask, [], id);
%Normalize with reference
referenceSpectrum = getReferenceSpectrum();
estimatedReflectanceColumn = reshape(estimatedReflectance, size(estimatedReflectance, 1), size(estimatedReflectance, 2)*size(estimatedReflectance, 3));
normalizedReflectance = estimatedReflectanceColumn ./ repmat(referenceSpectrum, 1, size(estimatedReflectanceColumn, 2));
normalizedReflectance = reshape(normalizedReflectance, size(estimatedReflectance, 1), size(estimatedReflectance, 2), size(estimatedReflectance, 3));
absorption = -log10(normalizedReflectance);
load(fullfile(getSetting('systemdir'), 'system.mat'), 'wavelength');
end

function newMap = removeInf(inMap)

origDim = size(inMap);
inMap = reshape(inMap, [origDim(1) * origDim(2), 1]);
maxMap = max(inMap(~isinf(inMap)));
maxLim = round(maxMap, -floor(log10(abs(maxMap))));

inMap(isinf(inMap)) = maxLim;

minMap = min(inMap(~(isinf(abs(inMap)) & (inMap < 0))));
minLim = round(minMap, -ceil(log10(abs(minMap))));

inMap((isinf(abs(inMap)) & (inMap < 0))) = minLim;

newMap = reshape(inMap, [origDim(1), origDim(2)]);
end
