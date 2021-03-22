function [T, measuredSpectra, adjustedSpectra, alphaCoeff] = evaluateColorchart(fileConditions, allowRoiSelection, selectedPatches, option)
%EVALUATECOLORCHART returns measured curves end evaluations in comparison
%to expected colorchart curves
%
%   [T, measuredSpectra, adjustedSpectra] = evaluateColorchart(
%   fileConditions, allowRoiSelection, selectedPatches) returns values for
%   evaluation of the colorchart
%
%   Input parameters
%   -fileConditions: conditions of HSI file to read
%   -allowRoiSelection: whether manual selection of colorchart ROI is
%   allowed
%   -selectedPatches: indexes of available patches to be used
%   -option: options for normalization when reading the hsi
%
%   Output parameters
%   -T: table including GOF, NMSE, RMSE of the evalueation of measured with
%   expected spectra
%   -measuredSpectra: measured spectra from the HSI
%   -adjustedSpectra: measured spectra from the HSI after adjustment
%   -alphaCoeff: coefficient used for level adjustment

experiment = getSetting('experiment');
patchOrder = getExpectedValues('colorchartOrder', experiment);
[expectedSpectra, patchNames, expectedWavelengths] = getExpectedValues();

if nargin < 2
    allowRoiSelection = false;
end

if nargin < 3 || isempty(selectedPatches)
    selectedPatches = 1:length(patchNames);
end

if nargin < 4
    option = [];
end
filename = getFilename(fileConditions{:});
[spectralData, ~, wavelengths] = loadH5Data(filename, experiment);
dispImage = getDisplayImage(spectralData, 'rgb');
[colorMasks, chartMask] = getColorchartMasks(dispImage, allowRoiSelection, filename);
actualSpectralVals = readHSI(spectralData, {chartMask, colorMasks}, option);

[reorderedSpectralVals, lineNames] = reorderSpectra(actualSpectralVals, patchOrder, patchNames, wavelengths, expectedWavelengths);
% [reorderedSpectralValsRaw, ~] = reorderSpectra(actualSpectralVals, chartColorOrder, spectraColorOrder, wavelengths, expectedWavelengths);

measuredSpectra = reorderedSpectralVals(selectedPatches, :);
standardSpectra = expectedSpectra(selectedPatches, :);

if ~strcmp(option, 'raw')
    [T, adjustedSpectra, alphaCoeff] = compareSpectra(standardSpectra, measuredSpectra, lineNames);
else
    plots(1, @plotColorChartSpectra, measuredSpectra, lineNames, 'measured-raw', [0, 0.005], false);
end

end