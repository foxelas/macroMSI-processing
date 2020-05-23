function [map] = getMap(msi, mapType)
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
% reference = getReference(getSetting('systemdir'), height, width);
reference = getReferenceFromMacbeth(height, width);
reference = raw2msi(reference, msiType);

switch mapType
    case 'opticalDensityMelanin'
        opticalDensity = double(log10(msi./reference));
        map = squeeze(opticalDensity(7, :, :));
        plotFunWrapper(fig, @plotMap, map, mapType);

    case 'opticalDensityHemoglobin'
        opticalDensity = double(log10(msi./reference));
        map = squeeze(opticalDensity(5, :, :)) - 1.15 .* squeeze(opticalDensity(7, :, :));
        plotFunWrapper(fig, @plotMap, map, mapType);

    case 'deepMelanin'
        normalizedReflectance = msi ./ reference;
        
    case 'totalMelanin'
        savedir = getSetting('savedir');
        mapdir = getSetting('map');
        setSetting('plotName', fullfile(savedir, mapdir, 'msi.png'));
        plotFunWrapper(2, @plotMSI, msi);
        normalizedReflectance = msi ./ reference;
        setSetting('plotName', fullfile(savedir, mapdir, 'normalizedMsi.png'));
        plotFunWrapper(3, @plotMSI, normalizedReflectance);
        absorption = -log10(normalizedReflectance);
        setSetting('plotName', fullfile(savedir, mapdir, 'absoprtion.png'));
        plotFunWrapper(4, @plotMSI, absorption);
        map = squeeze((absorption(5,:,:) - absorption(3,:,:)) ./ (575 - 505));
        setSetting('plotName', fullfile(savedir, mapdir, 'map.png'));
        figure(5);imshow(map); savePlot(5);
        setSetting('plotName', fullfile(savedir, mapdir, 'scaledMap.png'));
        figure(6); plotFunWrapper(6, @plotMap, map, [], false, 'Absorbance Slope (a.u.)');
    otherwise
end


end

% function slope = estimateSlope(x, y)
%         b = X\y
%         yCalc2 = X*b;
%         plot(x,yCalc2,'--')
%         legend('Data','Slope','Slope & Intercept','Location','best');
%         Rsq1 = 1 - sum((y - yCalc1).^2)/sum((y - mean(y)).^2)
% end 