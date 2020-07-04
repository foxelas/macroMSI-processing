function [melMap, hbMap] = getMap(raw, mapType, mask, id)
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

msiType = 'adjusted';
msi = raw2msi(raw, msiType);
[channels, height, width] = size(msi);

if nargin < 3
    mask = ones(height, width);
end 
% reference = getReference(getSetting('systemdir'), height, width);
reference = getReferenceFromMacbeth(height, width);
reference = raw2msi(reference, msiType);

savedir = getSetting('savedir');
mapdir = getSetting('map');
% setSetting('plotName', fullfile(savedir, mapdir, 'msi.png'));
% plotFunWrapper(2, @plotMSI, msi);
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
        logIm = log10(normalizedReflectance);
        opticalDensity = logIm;
        melMap = squeeze(opticalDensity(7, :, :));
        melBarTitle = 'Optical Density of Melanin';
        melSaveName = 'DingMel';
        melMapLimits = [-2.0, 0.2];
        hbMap = squeeze(opticalDensity(5, :, :)) - 1.15 .* squeeze(opticalDensity(7, :, :));
        hbBarTitle = 'Optical Density of Hemoglobin';
        hbMapSaveName = 'DingHb';
        hbMapLimits = [-3.0,1.0];

    case 'vasefi'
        logIm = log10(normalizedReflectance);
        absorption = -logIm; 
        fc = [450,465,505,525,575,605,630]';
        range = [605, 630];
        melMap = estimateSlope(absorption, range, fc);
        setSetting('plotName', fullfile(savedir, mapdir, 'scaledMap.png'));
        melBarTitle =  'Melanin Absorbance Slope (a.u.)';
        melSaveName = 'VasefiMel';
        melMapLimits = [-0.1, 0.02];
        range = [505, 575];
        hbMap = estimateSlope(absorption, range, fc);
        setSetting('plotName', fullfile(savedir, mapdir, 'scaledMap.png'));
        hbBarTitle =  'Hemoglobin Absorbance Slope (a.u.)';
        hbMapSaveName = 'VasefiHb';
        hbMapLimits = [-0.04,0.04];
    case 'diebele'
        logIm = log10(normalizedReflectance);
        opticalDensity = logIm;
        melMap = squeeze(opticalDensity(7, :, :));
        melBarTitle = 'Melanin Index';
        melSaveName = 'DiebeleMel';
        melMapLimits = [-2.0, 0.2];
        hbMap = squeeze(opticalDensity(7, :, :)) - squeeze(opticalDensity(5, :, :));
        hbBarTitle = 'Erythema Index';
        hbMapSaveName = 'DiebeleHb';
        hbMapLimits = [-3.0,1.0];
    case 'kapsokalyvas'
        %totalRed = squeeze(sum(squeeze(raw(:,:,:,3)), 1));
        %totalGreen = squeeze(sum(squeeze(raw(:,:,:,2)), 1));
        totalRed = squeeze(msi(7, :, :));
        totalGreen =  squeeze(msi(5, :, :));
        totalBlue = squeeze(msi(1, :, :));
        melMap = (totalRed - min(totalRed(:))) ./ (max(totalRed(:)) - min(totalRed(:)));
        melBarTitle = 'Melanin Contrast';
        melMap = (totalBlue - totalGreen) ./ (totalBlue + totalGreen);
        melBarTitle = 'Superficial Melanin';
        melSaveName = 'KapsokalyvasMel';
        melMapLimits = [0, 1];
        hbMap = (totalGreen - totalRed) ./ (totalGreen + totalRed);
        hbBarTitle = 'Hemoglobin Contrast';
        hbMapSaveName = 'KapsokalyvasHb';
        hbMapLimits = [-1, 1];
    case 'ours'
        [absorption, wavelength] = getAbsorptionByEstimatedReflectance(msi, mask, id);
        range = [605, 630];
        load(fullfile(getSetting('systemdir'), 'system.mat'), 'wavelength');
        melMap = estimateSlope(absorption, range, wavelength);
        melBarTitle =  'Melanin Absorbance Slope (a.u.)';
        melSaveName = 'oursMel';
        melMapLimits =  [];%[-1.505 * 0.001, -1.49 * 0.001];
        range = [505, 575];
        hbMap = estimateSlope(absorption, range, wavelength);
        hbBarTitle =  'Hemobglobin Absorbance Slope (a.u.)';
        hbMapSaveName = 'oursHb';
        hbMapLimits = []; %[1.45 * 0.001, 1.5 * 0.001];
        
    otherwise
        disp('Unsupported type')
end
setSetting('plotName', fullfile(savedir, mapdir, strcat(melSaveName,'ScaledMap.png')));
plotFunWrapper(1, @plotMap, melMap, mask, [], false, melBarTitle, melMapLimits);
setSetting('plotName', fullfile(savedir, mapdir, strcat(hbMapSaveName,'ScaledMap.png')));
plotFunWrapper(2, @plotMap, hbMap, mask, [], false, hbBarTitle,hbMapLimits);

end

function slopes = estimateSlope(image, range, fc)
        
rangeStart = range(1);
rangeEnd = range(2);
[b, m, n] = size(image);
condition = fc >= rangeStart & fc <= rangeEnd;
x = fc(condition);
columnImage = reshape(image, b, m * n);
y = columnImage(condition, :);
X = [ones(length(x),1) x];
b = X\y;
slopes = b(2,:);
slopes = reshape(slopes, m, n);

end 

function [absorption, wavelength] = getAbsorptionByEstimatedReflectance(I, mask, id)
setSetting('pixelValueSelectionMethod', 'adjusted');
setSetting('smoothingMatrixMethod', 'Cor_All');
setSetting('noiseType','fromOlympus');
estimatedReflectance = estimateReflectance(I, mask, [], id  );
%Normalize with reference 
referenceSpectrum = getReferenceSpectrum(); 
estimatedReflectanceColumn = reshape(estimatedReflectance, size(estimatedReflectance, 1),  size(estimatedReflectance, 2) *  size(estimatedReflectance, 3) );
normalizedReflectance = estimatedReflectanceColumn ./ repmat(referenceSpectrum, 1, size(estimatedReflectanceColumn, 2));
normalizedReflectance = reshape(normalizedReflectance, size(estimatedReflectance, 1),  size(estimatedReflectance, 2), size(estimatedReflectance, 3));
absorption = -log10(normalizedReflectance);
load(fullfile(getSetting('systemdir'), 'system.mat'), 'wavelength');
end
