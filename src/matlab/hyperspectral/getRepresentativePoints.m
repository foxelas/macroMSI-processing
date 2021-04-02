function [measured, curveNames] = getRepresentativePoints(input, xPoints, yPoints, limits, windowDim, lvlAdjCoeff)
%%GETREPRESENTATIVEPOINTS returns spectra of random points on the captured
%%surface.
% 
%     Arguments: 
%     input values for where to search for filename / value array
%     xPoints x coordinates for sample points 
%     yPoints y coordinates for sample points 
%     limits the limits for output spectrum plot 
%     windowDim window dimension if smooting is used 
%     lvlAdjCoeff multiplies for when level adjustment is used
% 
%     Usage: 
%     [measured, curveNames] = getRepresentativePoints(fileConditions,
%     xPoints, yPoints, limits, windowDim, lvlAdjCoeff)

if iscell(input)
    fileConditions = input;
    filename = getFilename(fileConditions{:});
    [raw, ~, ~] = loadH5Data(filename, getSetting('experiment'));
    targetName = fileConditions{4};
    spectralData = normalizeHSI(raw, targetName);
else 
    spectralData = input; 
    plotName = getSetting('plotName');
    splits = strsplit(plotName, '\');
    targetName = splits{end};
end 

if nargin < 2
    xPoints = 50:200:1088;
    yPoints = 50:200:982;
end

if nargin < 4
    limits = [0, 0.005];
end

if nargin < 5 || isempty(windowDim)
    windowDim = 0;
end

hasAdjustingCoeff = ~(nargin < 6);

if hasAdjustingCoeff
    if ~isequal(size(spectralData), size(lvlAdjCoeff))
        [m, n, ~] = size(spectralData);
        cropMask = getCaptureROImask(m,n);
        lvlAdjCoeff = lvlAdjCoeff(any(cropMask, 2), any(cropMask, 1), :);
        warning('Crop the image value: lvlAdjCoeff');
    end
    spectralData = spectralData .* lvlAdjCoeff;
end 


%% plot Y image with various points
baseImage = getDisplayImage(spectralData, 'rgb');

savedir = fullfile(getSetting('savedir'),  getSetting('saveFolder'));
setSetting('plotName', fullfile(savedir, strcat(targetName, '-points')));
plots(1, @plotPointsOnImage, baseImage, xPoints, yPoints, true);

wavelengths = getWavelengths(401); %without padding use (size(spectralData,3));
hasMultiple = ndims(spectralData) > 3;

%% Plot spectra

xLen = numel(xPoints);
yLen = numel(yPoints);

curveNames = cell(xLen*yLen, 1);
z = 0;
for i = xPoints
    for j = yPoints
        z = z + 1;
        curveNames{z} = sprintf('at (%d,%d)', i, j);
    end
end

if hasMultiple
    vals = zeros(xLen*yLen, numel(wavelengths), size(spectralData, 1));
    spectVals = zeros(size(spectralData, 1), xLen, yLen, numel(wavelengths));
else
    vals = zeros(xLen*yLen, numel(wavelengths));
    spectVals = zeros(xLen, yLen, numel(wavelengths));
end

z = 0;
for i = xPoints
    z = z + 1;
    y = 0;
    for j = yPoints
        y = y + 1;
        if hasMultiple
            for k = 1:size(spectralData, 1)
                curSpectra = squeeze(spectralData(k, z, y, :));
                curSpectra = pads(curSpectra);
                vals(sub2ind([xLen, yLen], z, y), :, k) = curSpectra;
                spectVals(k, z, y, :) = curSpectra;
            end
            
        else
            if windowDim <= 1
                curSpectra = squeeze(spectralData(i, j, :));
            elseif mod(windowDim, 2) == 0 % 2x2, 4x4
                windowX = (i - floor(windowDim/2)):(i + floor(windowDim/2) - 1);
                windowY = (j - floor(windowDim/2)):(j + floor(windowDim/2) - 1);
                curSpectra = mean(reshape(spectralData(windowX, windowY, :), [windowDim * windowDim, numel(wavelengths)]), 1)';
            else % 3x3
                windowX = (i - floor(windowDim/2)):(i + floor(windowDim/2));
                windowY = (j - floor(windowDim/2)):(j + floor(windowDim/2));
                curSpectra = mean(reshape(spectralData(windowX, windowY, :), [windowDim * windowDim, numel(wavelengths)]), 1)';
            end
            curSpectra = pads(curSpectra);
            spectVals(z, y, :) = curSpectra;
            vals(sub2ind([xLen, yLen], z, y), :) = curSpectra; %sub2ind([xLen, yLen], i, j)
        end
    end
end

if windowDim == 0
    curCase = 'measured-raw';
else
    curCase = strcat(num2str(windowDim), 'x', num2str(windowDim), 'window');
end

if max(vals(:)) > 0.15
   limits = [0, 1.2];
end 

plots(2, @plotColorChartSpectra, vals, curveNames, strcat(curCase, '_', targetName), ...
    limits, hasMultiple);

measured = reshape(spectVals, [1, numel(xPoints), numel(yPoints), numel(wavelengths)]);

end