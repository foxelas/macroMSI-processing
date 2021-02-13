function [measured, curveNames] = getRepresentativePoints(imgName, xPoints, yPoints, limits)

    filename = getFilename(getSetting('configuration'), imgName);
    [spectralData, ~, ~] = loadH5Data(filename, getSetting('experiment'));
    
    if nargin < 2
    xPoints = 50:200:1088;
    yPoints = 50:200:982;
    end 
    
    if nargin < 4 
        limits = [0,0.008];
    end 
        xLen = numel(xPoints);
    yLen = numel(yPoints);
    
    curveNames = cell(xLen * yLen, 1); 
    z = 0;
    for i = xPoints 
        for j = yPoints 
            z = z + 1;
            curveNames{z}= sprintf('at (%d,%d)', i,j);
        end 
    end 
    
    %% plot Y image with various points     
    setSetting('plotName',mkNewDir(getSetting('saveFolder'), imgName));
    baseImage = spectralData(:,:,100);
    plots(2, @plotYFromHSI, baseImage, 'Various points on the image', xPoints, yPoints);
    
    
    measured = reshape(showMultiplePointSpectra(k, xPoints, yPoints, 0, wavelengths', spectralData, limits), ...
        [1, numel(xPoints), numel(yPoints), numel(wavelengths)]);

end 