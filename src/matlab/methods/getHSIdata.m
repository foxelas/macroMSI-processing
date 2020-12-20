function spectralData = getHSIdata(raw, mask)

normByPixel =  getSetting('normByPixel');
spectralData = raw(any(mask, 2), any(mask, 1), :);

if normByPixel
    m = matfile(getSetting('normFilename'), 'Writable', false);
    normData = m.fullReflectanceByPixel;
    normData = normData(any(mask, 2), any(mask, 1), :);
    [n, m , d] = size(spectralData);
    spectralData = reshape(spectralData, [n * m, d]);
    fullReflectanceBypixel = reshape(normData, [n * m, d]);    
    spectralData = spectralData ./ fullReflectanceBypixel;
    spectralData(isnan(spectralData)) = 0;
    spectralData(isinf(spectralData)) = 0;
    spectralData = reshape(spectralData, [n, m, d]);
end 

end