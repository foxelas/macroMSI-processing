function [map] = getMap(msi, mapType)
%     GETMAP returns an optical map showing chromophore components of the macropathology image 
% 
%     Usage: 
%     [map] = getMap(msi)
%     [map] = getMap(msi, 'opticalDensityMelanin')

if varargin > 1
    mapType = 'opticalDensityMelanin';
end

[height, width, channels] = size(msi);
reference = getReference(getSetting('systemdir'), height, width);


switch mapType
    case 'opticalDensityMelanin'
        opticalDensity = double(log10(msi./reference));
        od630 = squeeze(opticalDensity(7, :, :));
        plotMap(od630, mapType);

    case 'opticalDensityHemoglobin'
        opticalDensity = double(log10(msi./reference));
        odhg = squeeze(opticalDensity(5, :, :)) - 1.15 .* squeeze(opticalDensity(7, :, :));
        plotMap(odhg, mapType);

    case 'deepMelanin'
        referenceReflectance =;
        normalizedReflectance = msi ./ reference;

    case 'totalMelanin'
    otherwise
end


end
