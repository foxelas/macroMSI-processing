
%% main
close all;
%clc;

%% PARAMETERS SET
%Modify userSettings.csv for options
userSettingsFile = '..\..\..\conf\hsiUserSettings.csv';
readWhite = false;

%% NEXT
setOpt(userSettingsFile);

dataDate = '20201111';
indir = strcat('F:\temp\mspi\saitama', dataDate, '_test/h5/');
setSetting('datadir', indir);
datafiles = dir(fullfile(indir, '*.h5'));

wavelengths = [380:780]';
m  = matfile('fullreflectance.mat', 'Writable', false);
if (~any(strcmp({whos(m).name}, 'maskWhite')) || readWhite)
    %For 99% reflectance
    [spectralData, imageXYZ, wavelengths] = loadH5Data(datafiles(1).name);
    figure(1);
    imagesc(imageXYZ(:, :, 2));
    title('Select ROI for Normalization Spectrum');
    maskWhite = roipoly;
    fullReflectance = getPatchValues(spectralData, maskWhite);
    save('fullreflectance.mat', 'fullReflectance', 'maskWhite');
end
load('fullreflectance.mat', 'fullReflectance', 'maskWhite');

setSetting('plotName', fullfile(getSetting('savedir'), 'whitePlot'));
plotFunWrapper(2, @plotSpectra, fullReflectance, wavelengths, '99%-white', 'Reflectance Spectrum of White Balance Sheet');

chartColorOrder = GetOrder();
expectedLabVals = GetExpectedLabValues();
expectedRGBVals = GetExpectedRGBValues();
[expectedSpectra, expectedWavelengths, spectraColorOrder] = GetExpectedSpectra();

function out = GetOrder()
outstruct = delimread('orderOfCurrentChart.txt', '\t', 'text');
out = outstruct.text;
end

function out = GetExpectedLabValues()
outstruct = delimread(fullfile(getSetting('systemdir'), 'ColorCheckerMicro_Matte_Lab_values.txt'), '\t', 'num');
out = outstruct.num;
end

function out = GetExpectedRGBValues()
outstruct = delimread(fullfile(getSetting('systemdir'), 'ColorCheckerMicro_Matte_RGB_values.txt'), '\t', 'num');
out = outstruct.num;
end

function [spectra, wavelengths, colorNames] = GetExpectedSpectra()
outstruct = delimread(fullfile(getSetting('systemdir'), 'ColorChecker_spectra.txt'), '\t', {'text', 'num'});
colorNames = outstruct.text;
colorNames = colorNames(2:length(colorNames));
wavelengths = outstruct.num(1, :);
spectra = outstruct.num(2:end, :);
end
