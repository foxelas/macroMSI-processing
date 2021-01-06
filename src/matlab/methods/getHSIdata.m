function spectralData = getHSIdata(raw, mask, settings)

if nargin < 3
    settings = 'useBlack';
end 

if ~strcmp(settings, 'raw')
    blackFilename = fullfile(getSetting('matdir'), getSetting('saveFolder'), 'blackReflectance.mat');
    load(blackFilename, 'blackReflectance');

    if (sum(size(blackReflectance) ~= size(raw)) >0 && (size(blackReflectance, 1) < size(raw,1)))
        [m,n,~] = size(raw);
        [x,y,~] = size(blackReflectance);
        raw = raw(floor(m/2)-floor(x/2):floor(m/2)+floor(x/2)-1, floor(n/2)-floor(y/2):floor(n/2)+floor(y/2)-1,:); 
        mask = mask(floor(m/2)-floor(x/2):floor(m/2)+floor(x/2)-1, floor(n/2)-floor(y/2):floor(n/2)+floor(y/2)-1);
    end 
end

spectralData = raw(any(mask, 2), any(mask, 1), :);
    
switch settings
    case 'useBlack'
        normByPixel =  getSetting('normByPixel');
  
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
        spectralData = (spectralData - blackReflectance) ./ normDenominator;
        
    case 'raw'
        spectralData = spectralData;
        
    case 'forExternalNormalization'               
        blackReflectance = blackReflectance(any(mask, 2), any(mask, 1), :);
        spectralData = (spectralData - blackReflectance);
        
    otherwise 
        error('Unsupported setting for normalization.');
end 

[m,n,z] = size(spectralData);
spectralData = reshape(spectralData, [m*n, z]);
spectralData(isnan(spectralData)) = 0;
spectralData(isinf(spectralData)) = 0;
spectralData = reshape(spectralData, [m,n,z]);

figure(6);imshow(squeeze(spectralData(:,:,200)));

end