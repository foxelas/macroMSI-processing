function spectralData = getHSIdata(raw, mask, settings, hasFilter)

if nargin < 3 || isempty(settings)
    settings = 'useBlack';
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

if ~strcmp(settings, 'raw')
    blackFilename = fullfile(getSetting('matdir'), configuration, strcat(num2str(integrationTime),  suffix,'_blackReflectance.mat')); %'basics.mat');
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
        
        wavelengths = size(spectralData,3);
        if ( wavelengths == 401)
            spectralData = (spectralData - blackReflectance) ./ normDenominator;
        elseif (wavelengths == 161)
            spectralData = (spectralData - blackReflectance(:,:,1:161)) ./ normDenominator(:,:,1:161);
        else
            spectralData = (spectralData - blackReflectance(:,:,162:end)) ./ normDenominator(:,:,162:end);
        end
        
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

figure(4);imshow(squeeze(spectralData(:,:,100)));

end