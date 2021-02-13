%% Compare effect of filter 
close all;

%% Settings 
userSettingsFile = '..\..\conf\hsiUserSettings.csv';

originDir = 'F:\temp\mspi';
dataDate = '20210111';
experiment = 'polarizing_filter_use_comparison';

%% Main 
setOpt(userSettingsFile);

matdir = fullfile(originDir, 'matfiles\hsi');
% Available configurations 'no_filter', 'filter'
configurations = {'no_filter', 'filter'}; 

[expected, patchNames] = getExpectedValues();

hasOldStruct = true; 
if ~hasOldStruct 
compData = struct('Experiment', [], 'Configuration', [], 'Position', [], 'Nmse', [], ...
    'Gof', [], 'Adjustment', [], 'MeanNmse', [], 'MeanGof' , [], 'NmseAdj', [], ...
    'GofAdj', [], 'MeanNmseAdj', [], 'MeanGofAdj', []);
end 

setSetting('saveImages', true);
for i = 1:2 
    configuration = configurations{i};
    patchOrder = getExpectedValues('colorchartOrder', configuration);

    indir = strcat(fullfile(originDir, '2_saitamaHSI\saitama'), dataDate, '_test\', configuration ,'\h5\');
    setSetting('datadir', indir);
    setSetting('saveFolder', fullfile(experiment, configuration));
    curSaveDir = mkNewDir(getSetting('savedir'), getSetting('saveFolder'));
    setSetting('matdir',  mkNewDir(matdir, getSetting('saveFolder')));

    datafiles = dir(fullfile(indir, '*.h5'));
    
    if hasOldStruct
         load(fullfile(curSaveDir, 'compData.mat'), 'compData');
    end 

%     blackFilename = datafiles( contains({datafiles.name}, 'black')).name;
%     [blackReflectance, imageXYZ, wavelengths] = loadH5Data(blackFilename, getSetting('saveFolder'));
%     %% plot Y image only
%     setSetting('plotName',mkNewDir(curSaveDir, strcat(configuration,'black')));
%     plots(1, @plotYFromHSI, imageXYZ, 'Luminance Y image');
%     
%     whiteFilename = datafiles( contains({datafiles.name}, 'white')).name;
%     [whiteReflectance, imageXYZ, ~] = loadH5Data(whiteFilename,  getSetting('saveFolder'));
%     
%     xPoints = 50:200:1088;
%     yPoints = 50:200:982;
%     %% plot Y image with various points       
%     setSetting('plotName',mkNewDir(curSaveDir, strcat(configuration, 'white')));
%     plots(2, @plotYFromHSI, imageXYZ, 'Various points on the image', xPoints, yPoints);
%     
%     whiteMinusBlackByPixel = bsxfun(@minus, whiteReflectance , blackReflectance);
%     save( fullfile(getSetting('matdir'), 'basics.mat'), 'blackReflectance', 'whiteReflectance', 'whiteMinusBlackByPixel', '-v7.3');   
% 
%     showMultiplePointSpectra(3, xPoints, yPoints, 1, wavelengths, whiteReflectance, configuration, curSaveDir);
%     
%     clear whiteReflectance; clear blackReflectance; 
%     
% %     calibrationFilename = datafiles( contains({datafiles.name}, 'calibration')).name;
% %     [spectralCalibration, imageXYZ, ~] = loadH5Data(calibrationFilename,  getSetting('saveFolder'));

    positions = {'left_up', 'left_down', 'right_down', 'right_up'};
    setSetting('normByPixel', true);
    for j = 1:4
        close all; 
        position = positions{j};
        
        positionFilename = datafiles( contains({datafiles.name}, position)).name;
        [spectralPosition, imageXYZ, ~] = loadH5Data(positionFilename,  getSetting('saveFolder'));
        
        %% plot Y image only
        setSetting('plotName',mkNewDir(curSaveDir, strcat(configuration,'calibration')));
        plots(3, @plotYFromHSI, imageXYZ, 'Luminance Y image');
        allowRoiSelection = true; 
        
        %% Temp because imageXYZ is corrupted 
        enableTemp = true;
        if enableTemp 
            for ik = 1:1376
            for jk = 1:1024
            I(ik,jk) = spectralPosition(ik,jk,200);
            end
            end
            imageXYZ = I; 
        end 
        %% End Temp 
        [colorMasks, chartMask] = getColorchartMasks(imageXYZ, allowRoiSelection, strcat(configuration, '_', num2str(j)));
        actualSpectralVals = readHSI(spectralPosition, {chartMask, colorMasks});
        [reorderedSpectralVals, lineNames] = reorderSpectra(actualSpectralVals, patchOrder, patchNames, wavelengths, expectedWavelengths);
        [reorderedSpectralValsRaw, ~] = reorderSpectra(actualSpectralVals, patchOrder, patchNames, wavelengths, expectedWavelengths);

        close all; 
       %% Standard (expected) color patch spectra
        plots(4, @plotColorChartSpectra, expected, lineNames, strcat('expected', '_', position));
        plots(5, @plotColorChartSpectra, reorderedSpectralVals, lineNames, strcat('measured', '_', position));
        
        differenceSpectra = expected - reorderedSpectralVals;
        plots(6, @plotColorChartSpectra, differenceSpectra, lineNames, strcat('difference', '_', position));

        gofs = applyRowFunc(reorderedSpectralVals, expected, @goodnessOfFit);
        nmses = applyRowFunc(reorderedSpectralVals, expected, @nmse);
        nn = 4 * (i-1) + j;
        compData(nn).Experiment = experiment;
        compData(nn).Configuration = configuration; 
        compData(nn).Position = position;
        compData(nn).Nmse = nmses; 
        compData(nn).Gof = gofs; 
        compData(nn).MeanNmse = mean(nmses); 
        compData(nn).MeanGof = mean(gofs);
        
       %% Standard (expected) color patch spectra after adjustment
        [adjustedReorderedSpectralVals, a] = adjustSpectra(reorderedSpectralVals, patchNames);

        plots(7, @plotColorChartSpectra, adjustedReorderedSpectralVals, lineNames, strcat('measured-adjusted', '_', position));

        differenceSpectraAdjusted = expected - adjustedReorderedSpectralVals;
        plots(8, @plotColorChartSpectra, differenceSpectraAdjusted, lineNames, strcat('difference-adjusted', '_', position));

        gofs2 = applyRowFunc(adjustedReorderedSpectralVals, expected, @goodnessOfFit);
        nmses2 = applyRowFunc(adjustedReorderedSpectralVals, expected, @nmse);

        compData(nn).NmseAdj = nmses2; 
        compData(nn).GofAdj = gofs2; 
        compData(nn).Adjustment = a; 
        compData(nn).MeanNmseAdj = mean(nmses2); 
        compData(nn).MeanGofAdj = mean(gofs2);      
    end 
    
    save(fullfile(curSaveDir, 'compData.mat'), 'compData');
    
end 

