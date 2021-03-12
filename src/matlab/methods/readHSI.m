function spectralData = readHSI(raw, mask, option, hasFilter)
%GETHSIDATA returns spectral data from HSI image
%
%   spectralData = readHSI(raw, mask, 'useBlack', false) returns a
%   cropped HSI
%
%   spectralData = readHSI(raw, {mainMask, subMasks}, 'useBlack', false)
%   returns the avarage for each submask
%

subMasks = [];
if iscell(mask)
    subMasks = mask{2:end};
    mask = mask{1};
end

if nargin < 3 || isempty(option)
    option = 'useBlack';
end

if nargin < 4
    hasFilter = true;
end

suffix = '';
if ~hasFilter
    suffix = 'NoFilter';
end

integrationTime = getSetting('integrationTime');
configuration = getSetting('configuration');

if ~strcmp(option, 'raw')
    blackFilename = fullfile(getSetting('matdir'), configuration, strcat(num2str(integrationTime), suffix, '_blackReflectance.mat')); %'basics.mat');
    load(blackFilename, 'blackReflectance');
    
    if (sum(size(blackReflectance) ~= size(raw)) > 0 && (size(blackReflectance, 1) < size(raw, 1)))
        [m, n, ~] = size(raw);
        [x, y, ~] = size(blackReflectance);
        raw = raw(floor(m/2)-floor(x/2):floor(m/2)+floor(x/2)-1, floor(n/2)-floor(y/2):floor(n/2)+floor(y/2)-1, :);
        mask = mask(floor(m/2)-floor(x/2):floor(m/2)+floor(x/2)-1, floor(n/2)-floor(y/2):floor(n/2)+floor(y/2)-1);
    end
end
%has opposite indexes because of any()
spectralData = raw(any(mask, 2), any(mask, 1), :);
clear 'raw';

switch option
    case 'useBlack'
        normByPixel = getSetting('normByPixel');
        
        if ~normByPixel
            load(blackFilename, 'whiteMinusBlack');
            normDenominator = whiteMinusBlack;
            clear whiteMinusBlack;
        else
            load(blackFilename, 'whiteMinusBlackByPixel');
            normDenominator = whiteMinusBlackByPixel;
            clear whiteMinusBlackByPixel;
        end
        
        blackReflectance = blackReflectance(any(mask, 2), any(mask, 1), :);
        normDenominator = normDenominator(any(mask, 2), any(mask, 1), :);
        
        wavelengths = size(spectralData, 3);
        if (wavelengths == 401)
            spectralData = (spectralData - blackReflectance) ./ normDenominator;
        elseif (wavelengths == 161)
            spectralData = (spectralData - blackReflectance(:, :, 1:161)) ./ normDenominator(:, :, 1:161);
        else
            spectralData = (spectralData - blackReflectance(:, :, 162:end)) ./ normDenominator(:, :, 162:end);
        end
        
    case 'raw'
        %do nothing
        
    case 'forExternalNormalization'
        blackReflectance = blackReflectance(any(mask, 2), any(mask, 1), :);
        spectralData = (spectralData - blackReflectance);
        
    otherwise
        error('Unsupported setting for normalization.');
end

[m, n, z] = size(spectralData);
spectralData = reshape(spectralData, [m * n, z]);
spectralData(isnan(spectralData)) = 0;
spectralData(isinf(spectralData)) = 0;
spectralData = reshape(spectralData, [m, n, z]);

% figure(4);imshow(squeeze(spectralData(:,:,100)));

if ~isempty(subMasks)
    y = size(subMasks, 3);
    target = spectralData;
    spectralData = zeros(y, z);
    
    for k = 1:y
        subMask = subMasks(:, :, k);
        %has opposite indexes because of any()
        maskedTarget = target(any(subMask, 2), any(subMask, 1), :);
        spectralData(k, :) = mean(reshape(maskedTarget, [size(maskedTarget, 1) * size(maskedTarget, 2), z]));
    end
    
end

end