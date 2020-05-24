function [map] = getMap(msi, mapType, mask, id)
%     GETMAP returns an optical map showing chromophore components of the macropathology image 
% 
%     Usage: 
%     [map] = getMap(msi)
%     [map] = getMap(msi, 'opticalDensityMelanin')

if nargin < 2
    mapType = 'opticalDensityMelanin';
end

msiType = 'adjusted';
msi = raw2msi(msi, msiType);
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
    case 'opticalDensityMelanin'
        map = squeeze(opticalDensity(7, :, :));
        barTitle = 'Optical Density of Melanin';
        saveName = 'DingMel';
        
    case 'opticalDensityHemoglobin'
        map = squeeze(opticalDensity(5, :, :)) - 1.15 .* squeeze(opticalDensity(7, :, :));
        barTitle = 'Optical Density of Hemoglobin';
        saveName = 'DingHb';

    case 'hemoglobin'
        fc = [450,465,505,525,575,605,630]';
        range = [505, 575];
        map = estimateSlope(absorption, range, fc);
        setSetting('plotName', fullfile(savedir, mapdir, 'scaledMap.png'));
        barTitle =  'Hemoglobin Absorbance Slope (a.u.)';
        saveName = 'VasefiHb';
        
    case 'deepMelanin'
        
    case 'totalMelanin'
        fc = [450,465,505,525,575,605,630]';
        range = [605, 630];
        map = estimateSlope(absorption, range, fc);
        setSetting('plotName', fullfile(savedir, mapdir, 'scaledMap.png'));
        barTitle =  'Melanin Absorbance Slope (a.u.)';
        saveName = 'VasefiMel';
    case 'mel'
        setSetting('pixelValueSelectionMethod', 'adjusted');
        setSetting('smoothingMatrixMethod', 'Cor_All');
        setSetting('noiseType','fromOlympus');
        estimatedReflectance = estimateReflectance(msi, mask, [], id  );
        absorption = -log10(estimatedReflectance);
        range = [605, 630];
        load(fullfile(getSetting('systemdir'), 'system.mat'), 'wavelength');
        map = estimateSlope(absorption, range, wavelength);
        barTitle =  'Melanin Absorbance Slope (a.u.)';
        saveName = 'ours';
    case 'hb'
        setSetting('pixelValueSelectionMethod', 'adjusted');
        setSetting('smoothingMatrixMethod', 'Cor_All');
        setSetting('noiseType','fromOlympus');
        estimatedReflectance = estimateReflectance(msi, mask, [], id  );
        absorption = -log10(estimatedReflectance);
        range = [505, 575];
        load(fullfile(getSetting('systemdir'), 'system.mat'), 'wavelength');
        map = estimateSlope(absorption, range, wavelength);
        barTitle =  'Hemobglobin Absorbance Slope (a.u.)';
        saveName = 'oursHb';
    otherwise
        disp('Unsupported type')
end
setSetting('plotName', fullfile(savedir, mapdir, strcat(saveName,'scaledMap.png')));
figure(6); plotFunWrapper(6, @plotMap, map, mask, [], false, barTitle);

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

