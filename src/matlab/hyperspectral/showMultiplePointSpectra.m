function spectVals = showMultiplePointSpectra(figNum, xPoints, yPoints, windowDim, wavelengths, spectralData, limits, suffix)
%SHOWMULTIPLEPOINTSPECTRA plots and returns pectra at specific points on 
%on the image
%
%   spectVals = showMultiplePointSpectra(figNum, xPoints, yPoints, 
%   windowDim, wavelengths, spectralData, limits, suffix)
%

    hasMultiple = ndims(spectralData) > 3;
    
    if nargin < 8
        suffix = '';
    end 
    
    %% plot Y image with various points       
    plotName = fullfile(getSetting('savedir'), getSetting('saveFolder'), strcat(figName, '-points', '.png'));
    setSetting('plotName', plotName);
    plots(figNum+1, @plotYFromHSI, imageXYZ, 'Various points on the image', xPoints, yPoints);

    %% Plot spectra
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
    
    if hasMultiple
        vals = zeros( xLen * yLen, numel(wavelengths), size(spectralData,1));
        spectVals = zeros(size(spectralData,1), xLen, yLen, numel(wavelengths));
    else 
        vals = zeros(xLen * yLen, numel(wavelengths));
        spectVals = zeros(xLen, yLen, numel(wavelengths));
    end 
    
    z = 0;
    for i = xPoints
        z = z+1;
        y = 0;
        for j = yPoints
            y = y+1;
            if hasMultiple 
                for k = 1:size(spectralData,1)
                    curSpectra =  squeeze(spectralData(k,z,y,:));
                    vals(sub2ind([xLen, yLen], z, y), :, k) = curSpectra;
                    spectVals(k,z,y,:) = curSpectra;
                end 
                
            else    
                if windowDim <= 1 
                    curSpectra =  squeeze(spectralData(i,j,:));
                elseif mod(windowDim,2) == 0 % 2x2, 4x4
                    windowX = (i-floor(windowDim/2)):(i+floor(windowDim/2)-1);
                    windowY = (j-floor(windowDim/2)):(j+floor(windowDim/2)-1);
                    curSpectra =  mean(reshape(spectralData(windowX,windowY,:), [windowDim * windowDim, numel(wavelengths)]), 1)';
                else % 3x3
                    windowX = (i-floor(windowDim/2)):(i+floor(windowDim/2));
                    windowY = (j-floor(windowDim/2)):(j+floor(windowDim/2));
                    curSpectra =  mean(reshape(spectralData(windowX,windowY,:), [windowDim * windowDim, numel(wavelengths)]), 1)';
                end
                spectVals(z, y, :) = curSpectra;
                vals(sub2ind([xLen, yLen], z, y), :) = curSpectra; %sub2ind([xLen, yLen], i, j)
            end
        end
    end

    if windowDim == 0 
        curCase = 'points';
    else 
        curCase = strcat(num2str(windowDim), 'x', num2str(windowDim), 'window');
    end 
    
    plots(figNum, @plotColorChartSpectra, vals, curveNames, strcat(curCase, '_', suffix), ...
        limits, hasMultiple);
end 