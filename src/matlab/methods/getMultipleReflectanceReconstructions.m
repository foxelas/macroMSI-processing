function [reconstructedArray, gfcArray, nmseArray] = getMultipleReflectanceReconstructions(msi, rgb, mask, measured, idd, methods)
% GETMULTIPLEREFLECTANCERECONSTRUCTIONS returns an array of reconstructed reflectances using different methods
%
% Usage:
% [reconstructedArray, gfcArray, nmseArray] = getMultipleReflectanceReconstructions(msi, rgb, mask, measured, idd, methods)

wavelength = 380:5:780;
lines = length(methods);
reconstructedArray = zeros(length(measured), lines);
gfcArray = zeros(1, lines);
nmseArray = zeros(1, lines);

for i = 1:lines
    method = methods{i};
    currentOptions = [];
    currentOptions.pixelValueSelectionMethod = getSetting('pixelValueSelectionMethod');
    currentOptions.noiseType = getSetting('noiseType');
    currentOptions.smoothingMatrixMethod = 'Cor_Sample';
    if strcmp(method, 'RGB-Simple')
        currentOptions.pixelValueSelectionMethod = 'rgb';
        currentOptions.noiseType = 'fromOlympus';
        img = readRGBImage(rgb);
        
        
    elseif strcmp(method, 'MSI-SpatioSpectral')
        currentOptions.rho = 0.6;
        currentOptions.windowDim = 3;
        currentOptions.noiseType = 'spatiospectralolympus';
        currentOptions.noiseParam = 0.0001;
        img = msi;
        
    elseif strcmp(method, 'MSI-Simple')
        currentOptions.noiseType = 'fromOlympus';
        currentOptions.noiseParam = 1;
        img = msi;
        
    elseif strcmp(method, 'MSI-Spatial')
        currentOptions.noiseType = 'spatial';
        currentOptions.noiseParam = [0.6310, 0.6310];
        img = msi;
    else
        error('Unavailable reconstructionMethod')
    end
    
    [reconstructedArray(:, i), gfcArray(i), nmseArray(i)] = estimateReflectance(img, mask, measured, idd, currentOptions);
    
end

if (getSetting('showImages'))
    spectra = [measured', reconstructedArray];
    lineNames = ['Measured', methods];
    setSetting('plotName', fullfile(getSetting('savedir'), getSetting('reflectanceEstimation'), getSetting('action'), ...
        strcat(getSetting('action'), '_', num2str(idd.Index))));
    plots(1, @plotReconstructedCurves, (spectra, lineNames, wavelength, 'Reflectance Estimation Comparison');
    pause(0.1)
end


end

function [rgb] = readRGBImage(whiteImg)
load('parameters\color_correction.mat', 'illuminant_gw1');
tempRGB = chromadapt(whiteImg, illuminant_gw1, 'ColorSpace', 'linear-rgb'); %color adjustment
rgb = reshape(tempRGB, [3, size(tempRGB, 1), size(tempRGB, 2)]);
end