function [outputArg1,outputArg2] = getMap(msi, mapType, saveOptions)

if varargin > 1 
    mapType = 'melanin';
end 

switch mapType
    case 'opticalDensityMelanin'
        opticalDensity = double(log10(msi ./ reference));
        od630 = squeeze(opticalDensity(7,:,:));
        saveOptions.figTitle = 'OD630nm - Melanin Map';
        saveOptions.relativeDir = 'od630';
        plotMap('od630', saveOptions);
        
    case 'opticalDensityHemoglobin'
        opticalDensity = double(log10(msi ./ reference));
        odhg = squeeze(opticalDensity(5,:,:)) - 1.15 .* squeeze(opticalDensity(7,:,:));
        saveOptions.figTitle = 'OD575nm - 1.15 OD630nm - Hemoglobin Map ';
        saveOptions.relativeDir = 'od630';
        plotMap('odhg', saveOptions);

    case 'deepMelanin'
    case 'totalMelanin'
    otherwise
end

        
end

