function spectralData = readHSI(raw, targetMask, option, hasFilter)
%GETHSIDATA returns spectral data from HSI image
%
%   spectralData = readHSI(raw, mask, 'useBlack', false) returns a
%   cropped HSI
%
%   spectralData = readHSI(raw, {mainMask, subMasks}, 'byPixel', false)
%   returns the avarage for each submask
%

subMasks = [];
if iscell(targetMask)
    subMasks = targetMask{2:end};
    targetMask = targetMask{1};
end

if nargin < 3 || isempty(option)
    option = getSetting('normalization');
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

[m, n, ~] = size(raw);
%for cases where size of image and size of white/black is different
mask = getCaptureROImask(m,n);

% targetMask = cropROI(targetMask, mask);
spectralData = cropROI(raw, targetMask);
[m, n, w] = size(spectralData);
clear 'raw';

whiteFilename = fullfile(getSetting('matdir'), configuration, strcat(num2str(integrationTime), suffix, '_fullReflectance.mat')); 
        
useBlack = true;
if useBlack && ~strcmp(option, 'raw')
    blackFilename = fullfile(getSetting('matdir'), configuration, strcat(num2str(integrationTime), suffix, '_blackReflectance.mat')); 
    load(blackFilename, 'blackReflectance');
    blackReflectance =  cropROI(blackReflectance, targetMask); 
end 

switch option        
    case 'raw'
        %do nothing
        useBlack = false;
        
    case 'byPixel'
        load(whiteFilename, 'fullReflectanceByPixel');
        whiteReflectance = cropROI(fullReflectanceByPixel, targetMask); 
        clear 'fullReflectanceByPixel';
        
    case 'uniSpectrum'
        load(whiteFilename, 'uniSpectrum');
        whiteReflectance = reshape(repmat(uniSpectrum, m*n, 1), m, n, w);
        
    case 'bandmax'
        bandMaxSpectrum = max(reshape(spectralData, m * n, w), [], 1);
        setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('experiment'), getSetting('saveFolder'), 'bandmax'));
        plots(1, @plotSpectra, bandMaxSpectrum, getWavelengths(w), 'Bandmax spectrum', 'Bandmax Spectrum for the current Image');
        whiteReflectance = reshape(repmat(bandMaxSpectrum, m*n, 1), m, n, w);
        
    case 'forExternalNormalization'
        useBlack = false;
        spectralData = (spectralData - blackReflectance);
        
    otherwise
        error('Unsupported setting for normalization.');
end

if useBlack
    normDenominator = whiteReflectance - blackReflectance;
    spectralData = (spectralData - blackReflectance) ./ normDenominator;
end   

spectralData = max(spectralData, 0);
spectralData(isnan(spectralData)) = 0;
spectralData(isinf(spectralData)) = 0;

% figure(4);imshow(squeeze(spectralData(:,:,100)));

if ~isempty(subMasks)
    y = size(subMasks, 3);
    target = spectralData;
    spectralData = zeros(y, w);
    
    isSinglePoint = sum(subMasks(:,:,1), 'all') == 1;
    for k = 1:y
        subMask = subMasks(:, :, k);
        patchSpectra = cropROI(target, subMask); 
        if isSinglePoint 
            spectralData(k,:) = patchSpectra;
        else 
            spectralData(k, :) = mean(reshape(patchSpectra, [size(patchSpectra, 1) * size(patchSpectra, 2), w]));
        end 
    end
else 
    target = spectralData;
    spectralData = zeros(1, w);
    spectralData(1, :) = mean(reshape(target, [size(target, 1) * size(target, 2), w]));  
end

end

function hsiOut = cropROI(hsiIn, cropMask)
    %has opposite indexes because of any()
    hsiOut = hsiIn(any(cropMask, 2), any(cropMask, 1), :);
end 