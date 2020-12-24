function spectralData = getHSIdata(raw, mask)

normByPixel =  getSetting('normByPixel');
useBlack = true; 

if useBlack 
    blackFilename = fullfile(getSetting('matdir'), getSetting('saveFolder'), 'blackReflectance.mat');
    load(blackFilename, 'blackReflectance');
    if ~normByPixel 
        load(blackFilename, 'whiteMinusBlack');
        normDenominator = whiteMinusBlack;
        clear whiteMinusBlack;
    else 
        load(blackFilename, 'whiteMinusBlackByPixel');
        normDenominator = whiteMinusBlackByPixel;
        clear whiteMinusBlackByPixel; 
    end 
        
    if (sum(size(blackReflectance) ~= size(raw)) >0 && (size(blackReflectance, 1) < size(raw,1)))
        [m,n,~] = size(raw);
        [x,y,~] = size(blackReflectance);
        raw = raw(floor(m/2)-floor(x/2):floor(m/2)+floor(x/2)-1, floor(n/2)-floor(y/2):floor(n/2)+floor(y/2)-1,:); 
        mask = mask(floor(m/2)-floor(x/2):floor(m/2)+floor(x/2)-1, floor(n/2)-floor(y/2):floor(n/2)+floor(y/2)-1);
    end 
end

spectralData = raw(any(mask, 2), any(mask, 1), :);

if useBlack 
    blackReflectance = blackReflectance(any(mask, 2), any(mask, 1), :);
    normDenominator = normDenominator(any(mask, 2), any(mask, 1), :);
    spectralData = (spectralData - blackReflectance) ./ normDenominator;
end

[m,n,z] = size(spectralData);
spectralData = reshape(spectralData, [m*n, z]);
spectralData(isnan(spectralData)) = 0;
spectralData(isinf(spectralData)) = 0;
spectralData = reshape(spectralData, [m,n,z]);

figure(6);imshow(squeeze(spectralData(:,:,200)));

end