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

if isstruct(id)
    id = id.Index;
end

% reference = getReference(getSetting('systemdir'), height, width);
reference = getReferenceFromMacbeth(height, width);
reference = raw2msi(reference, msiType);

savedir = getSetting('savedir');
mapdir = getSetting('map');
% setSetting('plotName', fullfile(savedir, mapdir, 'msi.png'));
% plotFunWrapper(2, @plotMSI, msi);
msi = reshape(msi, [channels * height * width, 1]);
msi(msi == 0) =  0.0000000001;
msi = reshape(msi, [channels, height, width]);
mask3d = permute(repmat(mask, [1, 1, size(msi, 1)]), [3, 1, 2]);
msi(~mask3d) = nan;
nmsi = msi ./ reference;
% setSetting('plotName', fullfile(savedir, mapdir, 'normalizedMsi.png'));
% plotFunWrapper(3, @plotMSI, normalizedReflectance);

logIm = log10(nmsi);
opticalDensity = logIm;
absorption = -logIm;
% setSetting('plotName', fullfile(savedir, mapdir, 'absoprtion.png'));
% plotFunWrapper(4, @plotMSI, absorption);
switch lower(mapType)
    case 'ding'
        melMap = squeeze(opticalDensity(7, :, :));
        melBarTitle = 'Optical Density of Melanin (a.u.)';
        melSaveName = 'DingMel';
        melMapLimits = [-1.5, 0.5];
        hbMap = squeeze(opticalDensity(5, :, :)) - 1.15 .* squeeze(opticalDensity(7, :, :));
        hbBarTitle = 'Optical Density of Hemoglobin (a.u.)';
        hbMapSaveName = 'DingHb';
        hbMapLimits = [-5, 1];

    case 'vasefi'
        fc = [450, 465, 505, 525, 575, 605, 630]';
        melMap = getSlope(absorption, find(fc == 605), find(fc == 630));
        setSetting('plotName', fullfile(savedir, mapdir, 'scaledMap.png'));
        melBarTitle = 'Relative Melanin Concentration (a.u.)';
        melSaveName = 'VasefiMel';
        melMapLimits = [];
        
        [cHbO, cHbR] = estimageLR(absorption, find(fc == 505), find(fc == 575));
        hbMap = cHbO + cHbR;
        setSetting('plotName', fullfile(savedir, mapdir, 'scaledMap.png'));
        hbBarTitle = 'Relative Total Hemoglobin Concentration (a.u.)';
        hbMapSaveName = 'VasefiHb';
        hbMapLimits = [];
        
        cHbO = mat2gray(removeInf(cHbO));
        setSetting('plotName', fullfile(savedir, mapdir, strcat('VasefiHbO', '_', num2str(id), '_', 'ScaledMap.png')));
        plotFunWrapper(3, @plotMap, cHbO, mask, [], false, 'Relative HbO Concentration (a.u.)', []);
        cHbR = mat2gray(removeInf(cHbR));
        setSetting('plotName', fullfile(savedir, mapdir, strcat('VasefiHbR', '_', num2str(id), '_', 'ScaledMap.png')));
        plotFunWrapper(4, @plotMap, cHbR, mask, [], false, 'Relative HbR Concentration (a.u.)', []);
        

    case 'diebele'
        melMap = 100 * (squeeze(absorption(6, :, :)) - squeeze(absorption(7, :, :)));
        melBarTitle = 'Melanin Index (a.u.)';
        melSaveName = 'DiebeleMel';
        melMapLimits = [];
        hbMap = 100 * (squeeze(absorption(4, :, :)) - squeeze(absorption(7, :, :)));
        hbBarTitle = 'Erythema Index (a.u.)';
        hbMapSaveName = 'DiebeleHb';
        hbMapLimits = [];

    case 'kapsokalyvas'
        %totalRed = squeeze(sum(squeeze(raw(:,:,:,3)), 1));
        %totalGreen = squeeze(sum(squeeze(raw(:,:,:,2)), 1));
        totalRed = squeeze(nmsi(7, :, :));
        totalGreen = squeeze(nmsi(5, :, :));
        totalBlue = squeeze(nmsi(1, :, :));
        melMap = (-1) * (totalRed - min(totalRed(:))) ./ (max(totalRed(:)) - min(totalRed(:)));
        melBarTitle = 'Melanin Homogeneity (a.u.)';
        melSaveName = 'KapsokalyvasMel';
        melMapLimits = [];
        
        supMelMap = (-1) * (totalBlue - totalGreen) ./ (totalBlue + totalGreen);
        supMelMap = mat2gray(removeInf(supMelMap));
        setSetting('plotName', fullfile(savedir, mapdir, strcat('KapsokalyvasSupMel', '_', num2str(id), '_', 'ScaledMap.png')));
        plotFunWrapper(3, @plotMap, supMelMap, mask, [], false, 'Superficial Melanin Homogeneity (a.u.)', []);
            
        hbMap = (-1) * (totalGreen - totalRed) ./ (totalGreen + totalRed);
        hbBarTitle = 'Hemoglobin Homogeneity (a.u.)';
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
        melBarTitle = 'Melanin Index (a.u.)';
        melSaveName = 'KuzminaMel';
        melMapLimits = [0, 25];
        calibCoeff = redCoeff / greenCoeff;
        hbMap = calibCoeff * totalRed ./ totalGreen;
        hbBarTitle = 'Hemoglobin Index (a.u.)';
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

melMap = mat2gray(removeInf(melMap));
hbMap = mat2gray(removeInf(hbMap));
setSetting('plotName', fullfile(savedir, mapdir, strcat(melSaveName, '_', num2str(id), '_', 'ScaledMap.png')));
plotFunWrapper(1, @plotMap, melMap, mask, [], false, melBarTitle, melMapLimits);
setSetting('plotName', fullfile(savedir, mapdir, strcat(hbMapSaveName, '_', num2str(id), '_', 'ScaledMap.png')));
plotFunWrapper(2, @plotMap, hbMap, mask, [], false, hbBarTitle, hbMapLimits);

end

function slope = getSlope(image, l1, l2)

slope = (squeeze(image(l2,:,:)) - squeeze(image(l1,:,:))) / (l2 - l1);

end 

function [cHbO, cHbR] = estimageLR(image, l1, l2, fc)

[b, m, n] = size(image);
if nargin < 4
    condition = zeros(b,1);
    condition(l1:l2) = 1; 
    condition = logical(condition);
    x1 = [19946;32496.4000000000;55540];
    x2 = [23774.4000000000;35944;40092]; 
else
    rangeStart = l1;
    rangeEnd = l2;
    condition = fc >= rangeStart & fc <= rangeEnd;
    load('parameters\extinctionCoefficients.mat','extCoeffHbO', 'extCoeffHbR', 'hbLambda');
    rng = fc(condition);
    rngIdxs = arrayfun(@(z) find((hbLambda - z) >= 0 , 1), rng);
    x1 = extCoeffHbO(rngIdxs);
    x2 = extCoeffHbR(rngIdxs);
   
end

%X = [ones(length(x), 1), x];
X = [x1, x2]; 
columnImage = reshape(image, b, m*n);
y = columnImage(condition, :);

b = X \ y;
cHbO = b(1,:);
cHbR = b(2, :);
cHbO = reshape(cHbO, m, n);
cHbR = reshape(cHbR, m, n);


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
rounding = -floor(log10(abs(maxMap)));
if isinf(rounding)
    rounding = 0; 
end 
maxLim = round(maxMap, rounding);

inMap(isinf(inMap)) = maxLim;

minMap = min(inMap(~(isinf(abs(inMap)) & (inMap < 0))));
rounding = --ceil(log10(abs(minMap)));
if isinf(rounding)
    rounding = 0; 
end 
minLim = round(minMap, rounding);

inMap((isinf(abs(inMap)) & (inMap < 0))) = minLim;

newMap = reshape(inMap, [origDim(1), origDim(2)]);
end
