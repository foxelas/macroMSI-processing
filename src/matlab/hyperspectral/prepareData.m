
%% main
close all;
clear all;
%clc;

%% PARAMETERS SET
%Modify userSettings.csv for options
userSettingsFile = '..\..\conf\hsiUserSettings.csv';
readWhite = true;
readBlack = true;
wavelengths = [380:780]';

% Available dates,  '20201209', '20201111', '20201218', '20201220'
dataDate = '20201220';
% Available configurations 'singleLightFar', 'singleLightClose',
% 'doubleLightClose', 'noLight'
configuration = 'doubleLightClose'; 

%% NEXT
setOpt(userSettingsFile);

indir = strcat('F:\temp\mspi\saitama', dataDate, '_test/h5/');
setSetting('datadir', indir);
% datafiles = dir(fullfile(indir, '*.h5'));

matfileDir = getSetting('matdir');

%% Settings for normalization 1. Single Value, 2. PixelByPixel
setSetting('normByPixel', true);
normByPixel = getSetting('normByPixel');
hasSmoothing = false;

if ~normByPixel 
    setSetting('normFilename', mkNewDir(matfileDir, configuration, 'fullreflectance.mat')); 
    setSetting('saveFolder', strcat(configuration, '_byAverage'));
else 
    setSetting('normFilename', mkNewDir(matfileDir, configuration, 'fullreflectance_ByPixel.mat'));
    if hasSmoothing
        setSetting('saveFolder', strcat(configuration, '_byPixel_withSmoothing'));
    else 
        setSetting('saveFolder',strcat(configuration, '_byPixel_noSmoothing'));
    end 
end 
normFilename = getSetting('normFilename');
saveFolder = getSetting('saveFolder');

%% Values for normalization
if readWhite
    filename = getFilename(configuration, 'whiteReflectance');
    %For 99% reflectance
    [spectralData, imageXYZ, wavelengths] = loadH5Data(filename, configuration);
    whiteSize = size(spectralData);
    figure(1);
    dispImage = imageXYZ(:, :, 2);
    imagesc(dispImage);
    if ~normByPixel 
        title('Select ROI for Normalization Spectrum');
        maskWhite = roipoly;
        fullReflectance = getPatchValues(spectralData, maskWhite);
        save(normFilename, 'fullReflectance', 'maskWhite');
        setSetting('plotName', mkNewDir(getSetting('savedir'), saveFolder, 'whitePlot_'));
        plotFunWrapper(2, @plotSpectra, fullReflectance, wavelengths, '99%-white', 'Reflectance Spectrum of White Balance Sheet');
    else 
        xPoints = [100, 500, 800]; 
        if strcmp(configuration, 'singleLightFar')
            xPoints = xPoints + 500;
        end 
        if strcmp(configuration, 'singleLightClose')
            xPoints = xPoints + 200;
        end 
        yPoints = [100, 500, 800] ;
        fullReflectanceByPixel = spectralData; 
        pointSpectra = reshape(fullReflectanceByPixel(xPoints,yPoints,:), [3*3, length(wavelengths)]);
        setSetting('plotName', mkNewDir(getSetting('savedir'), saveFolder, 'whitePlot_by_pixel_unsmoothened'));
        plotFunWrapper(2, @plotSpectra, pointSpectra,...
            wavelengths, cellfun(@num2str, num2cell(1:9), 'UniformOutput', false), 'Reflectance Spectrum of various pixels in White Balance Sheet');
        
        fullReflectanceByPixel = smoothImage(hasSmoothing, fullReflectanceByPixel, wavelengths);
        pointSpectra = reshape(fullReflectanceByPixel(xPoints,yPoints,:), [3*3, length(wavelengths)]);
        setSetting('plotName',mkNewDir(getSetting('savedir'), getSetting('saveFolder'), 'whitePlot_by_pixel'));
        plotFunWrapper(2, @plotSpectra, pointSpectra,...
            wavelengths, cellfun(@num2str, num2cell(1:9), 'UniformOutput', false), 'Reflectance Spectrum of various pixels in White Balance Sheet');
        
        if sum(pointSpectra(2,:)' == squeeze(fullReflectanceByPixel(xPoints(2), yPoints(1), :))) ~= 401
            error('Not coinciding order of points');
        end 
        
        [ yy , xx] = meshgrid(xPoints,yPoints); 
        xx = xx(:);
        yy = yy(:);
        
        fig1 = figure(1);
        imagesc(dispImage);
        hold on;
        for i = 1:length(xx)
            plot(xx(i),yy(i),'rx', 'MarkerSize', 20, 'LineWidth', 5);
            textStr = sprintf('P%d at (%d,%d)', i,xx(i), yy(i));
            text(xx(i),yy(i), textStr);
        end 
        hold off;
        setSetting('plotName',mkNewDir(getSetting('savedir'), getSetting('saveFolder'), 'whitePlot_by_pixel-points'));
        savePlot(fig1);
        
        %save(normFilename, 'fullReflectanceByPixel', '-v7.3');
    end 
end

if readBlack 
    filename = getFilename('noLight', 'capOn');
    %For 99% reflectance
    [spectralData, imageXYZ, wavelengths] = loadH5Data(filename, configuration);
    figure(3);
    imagesc(imageXYZ(:, :, 2));
    blackReflectance = spectralData; 
    blackReflectance = smoothImage(hasSmoothing, blackReflectance, wavelengths);

    
    [x,y,z] = size(blackReflectance);
    if (sum(size(blackReflectance) ~= whiteSize) >0 && (size(blackReflectance, 1) < whiteSize(1)) && normByPixel)
        m = whiteSize(1);
        n = whiteSize(2);
        fullReflectanceByPixel = fullReflectanceByPixel(floor(m/2)-floor(x/2)+1:floor(m/2)+floor(x/2), floor(n/2)-floor(y/2)+1:floor(n/2)+floor(y/2),:); 
    end
        
    if ~normByPixel 
        whiteMinusBlack = fullReflectance - reshape(blackReflectance, [x*y, z]);
        whiteMinusBlack = reshape(whiteMinusBlack, [x,y,z]);
        save( mkNewDir(matfileDir, saveFolder, 'blackReflectance.mat'), 'blackReflectance', 'whiteMinusBlack');    
    else  
        whiteMinusBlackByPixel = fullReflectanceByPixel - blackReflectance;
        save( mkNewDir(matfileDir, saveFolder, 'blackReflectance.mat'), 'blackReflectance', 'whiteMinusBlackByPixel');    
    end 
end 

chartColorOrder = GetOrder(configuration);
expectedLabVals = GetExpectedLabValues();
expectedRGBVals = GetExpectedRGBValues();
[expectedSpectra, expectedWavelengths, spectraColorOrder] = GetExpectedSpectra();

function out = GetOrder(configuration)
outstruct = delimread(fullfile(getSetting('datasetSettingsDir'), strcat(configuration, 'PatchOrder.txt')), '\t', 'text');
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

function searchName = getFilename(configuration, content)
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    dataTable = readtable( fullfile(getSetting('datasetSettingsDir'), 'dataSetCharacteristics.csv')); 
    setId = strcmp(dataTable.configuration, configuration) & strcmp(dataTable.content, content);
    searchName = dataTable.filename{setId};
end 

function result = smoothImage(hasSmoothing, target, wavelengths)
    if hasSmoothing
        % smoothing the data with moving average filter 
        windowSize = 10; 
        b = (1/windowSize)*ones(1,windowSize);
        a = 1;
        y = filter(b,a,target, [], 1);
        result = filter(b,a,y, [], 2);
    else
        result = target;
    end
        
end