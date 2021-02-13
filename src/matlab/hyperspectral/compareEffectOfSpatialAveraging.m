%% Compare effect of spatial averaging 
close all;

%% Settings 
userSettingsFile = '..\..\conf\hsiUserSettings.csv';

originDir = 'F:\temp\mspi';
dataDate = '20201220';
% Available configurations 'singleLightFar', 'singleLightClose',
% 'doubleLightClose', 'noLight'
configuration = 'doubleLightClose'; 
experiment = 'spatial_average_comparison';

%% Main 
setOpt(userSettingsFile);

indir = strcat(fullfile(originDir, '2_saitamaHSI\'), dataDate, '_test\h5\');
setSetting('datadir', indir);
matdir = fullfile(originDir, 'matfiles\hsi');
setSetting('matdir', matdir);
setSetting('saveFolder', experiment);

curSaveDir = mkNewDir(getSetting('savedir'), getSetting('saveFolder'));

names = {'doubleLightClose', 'singleLightClose', 'singleLightFar'};
spectVals = zeros(numel(names), 6, 5, 401);
spect2Vals = zeros(numel(names), 6, 5, 401);
spect3Vals = zeros(numel(names), 6, 5, 401);

for k =1:numel(names)
    filename = getFilename(names{k}, 'colorchart');
    tempDataDate = strsplit(filename, '_');
    indir = strcat(fullfile(originDir, '2_saitamaHSI\'), tempDataDate{1}, '_test\h5\');
    setSetting('datadir', indir);
    [spectralData, imageXYZ, wavelengths] = loadH5Data(filename, experiment);
    
    xPoints = 50:200:1088;
    yPoints = 50:200:982;
    
    %% plot Y image only
    setSetting('plotName',mkNewDir(curSaveDir, strcat(names{k} ,'-white')));
    plots(1, @plotYFromHSI, imageXYZ, 'Luminance Y image');
    
    %% plot Y image with various points     
    setSetting('plotName',mkNewDir(curSaveDir, strcat(names{k} ,'-white')));
    plots(2, @plotYFromHSI, imageXYZ, 'Various points on the image', xPoints, yPoints);
    
    spectVals(k, :, :, :) = reshape(showMultiplePointSpectra(3, xPoints, yPoints, 1, wavelengths, spectralData, names{k}, curSaveDir), [1, numel(xPoints), numel(yPoints), numel(wavelengths)]);
    spect2Vals(k, :, :, :) = reshape(showMultiplePointSpectra(4, xPoints, yPoints, 2, wavelengths, spectralData, names{k}, curSaveDir), [1, numel(xPoints), numel(yPoints), numel(wavelengths)]);
    spect3Vals(k, :, :, :) = reshape(showMultiplePointSpectra(5, xPoints, yPoints, 3, wavelengths, spectralData, names{k}, curSaveDir), [1, numel(xPoints), numel(yPoints), numel(wavelengths)]);

end 
    
load('parameters/solax_reconstructed_spectrum.mat', 'solaxSpec', 'solaxLocalMaxWavelengths');
nn = numel(xPoints);
mm = numel(yPoints);
allCurves = zeros(3, 3, 401);
allCurves(3,1,:) = removePoints(spectVals(3,:,:,:), nn, mm, [1:10, 18,19, 23, 24]);
allCurves(3,2,:) = removePoints(spect2Vals(3,:,:,:), nn, mm, [1:10, 18,19, 23, 24]);
allCurves(3,3,:) = removePoints(spect3Vals(3,:,:,:), nn, mm, [1:10, 18,19, 23, 24]);
allCurves(2,1,:) = removePoints(spectVals(2,:,:,:), nn, mm, [13, 18]);
allCurves(2,2,:) = removePoints(spect2Vals(2,:,:,:), nn, mm, [13, 18]);
allCurves(2,3,:) = removePoints(spect3Vals(2,:,:,:), nn, mm, [13, 18]);
allCurves(1,1,:) = removePoints(spectVals(1,:,:,:), nn, mm, [13, 18]);
allCurves(1,2,:) = removePoints(spect2Vals(1,:,:,:), nn, mm, [13, 18]);
allCurves(1,3,:) = removePoints(spect3Vals(1,:,:,:), nn, mm, [13, 18]);

figure(5);clf;
hold on 
k = 0;
h = zeros(9,1);
for i = 1:3
    for j =1:3
        k = k + 1;
        if j == 1
            mark = '-';
        elseif j == 2
            mark = '--';
        else 
            mark = '.';
        end 
        dispName = sprintf('%s + %dx%d', names{i}, j, j);
        h(k) = plot(wavelengths, squeeze(allCurves(i,j,:)), mark, 'DisplayName', dispName, 'LineWidth', 3);
    end 
end 
h(10) = plot(wavelengths, solaxSpec ./ 10^(5), 'DisplayName', 'Solax-iO illumination spectrum (div by 10^5)', 'LineWidth', 5);
for i = 1:numel(solaxLocalMaxWavelengths)
    xline(solaxLocalMaxWavelengths(i), '-', sprintf('Solax-iO Local Max at %d nm', solaxLocalMaxWavelengths(i)), 'LineWidth', 3);
end 
hold off 
ylim([0,0.005]);
title('Spectra mean for selected spectra with spatial averaging');
legend(h, 'Location', 'eastoutside');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
setSetting('plotName',mkNewDir(curSaveDir, strcat(experiment ,'-pointSpectraMean')));
savePlot(5);

function [res] = removePoints(vals, nn, mm, toRemove)
valsTemp = reshape(squeeze(vals), [nn*mm, 401]);
valsTemp2 = zeros(nn*mm-numel(toRemove), 401);
j = 0;
for i = 1:size(valsTemp,1)
    if isempty(find(toRemove == i))
        j = j + 1;
        valsTemp2(j, :) = valsTemp(i,:);
    end 
end
res = squeeze(mean(valsTemp2, 1));
end 