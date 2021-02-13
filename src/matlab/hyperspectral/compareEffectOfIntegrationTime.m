%% Compare different exposure times
close all;

userSettingsFile = '..\..\conf\hsiUserSettings.csv';

dataDate = 'saitama20201209';
originDir = 'F:\temp\mspi';
exp100name = '20201209_173035-manualexposure100.h5';
exp200name = '20201209_171535-manualexposure200.h5';
experiment = 'integration_time_comparison';

setOpt(userSettingsFile);

indir = strcat(fullfile(originDir, '2_saitamaHSI\'), dataDate, '_test\h5\');
setSetting('datadir', indir);
matdir = fullfile(originDir, 'matfiles\hsi');
setSetting('matdir', matdir);
setSetting('saveFolder', experiment);

names = {exp100name, exp200name};
altNames = {'int100', 'int200'};
curSaveDir = mkNewDir(getSetting('savedir'), getSetting('saveFolder'));

spectVals = zeros(2, 7, 7, 401);
for k =1:2
    [spectralData, imageXYZ, wavelengths] = loadH5Data(names{k}, experiment);
    
    xPoints = 50:200:size(spectralData, 1);
    yPoints = 50:200:size(spectralData, 2);
    
    %% plot Y image only
    setSetting('plotName',mkNewDir(curSaveDir, strcat(altNames{k} ,'-white')));
    plots(1, @plotYFromHSI, imageXYZ, 'Luminance Y image');
    
    %% plot Y image with various points     
    setSetting('plotName',mkNewDir(curSaveDir, strcat(altNames{k} ,'-white')));
    plots(2, @plotYFromHSI, imageXYZ, 'Various points on the image', xPoints, yPoints);
    
    figure(3); clf; 
    z = 0;
    for i = xPoints
        z = z+1;
        y = 0;
        for j = yPoints
            y = y+1;
            spectVals(k, z, y, :) = squeeze(spectralData(i,j,:)); %sub2ind([numel(xPoints), numel(yPoints)], i, j)
            hold on
            plot(wavelengths, squeeze(spectralData(i,j,:)), 'DisplayName', sprintf('at (%d,%d)', i,j));
            hold off
        end
    end

    ylim([0,0.0012]);
    title('Spectra from random points across all the area of the white image ');
    legend('Location', 'eastoutside');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    setSetting('plotName',mkNewDir(curSaveDir, strcat(altNames{k} ,'-pointSpectra')));
    savePlot(3);
end 

figure(4); clf; 
for i = 1:numel(xPoints)
    for j = 1:numel(yPoints)
        diffSpect = squeeze(spectVals(1,i,j,:)) - squeeze(spectVals(2,i,j,:));
        hold on
        plot(wavelengths, diffSpect, 'DisplayName', sprintf('at (%d,%d)', i,j));
        hold off
    end
end 
hold on
plot(wavelengths, squeeze(mean(mean(spectVals(1,:,:,:) - spectVals(2,i,j,:), 2), 3)), 'DisplayName', 'Dif Mean for 100ms - 200ms', 'LineWidth', 5);
hold off
ylim([-0.0005, 0.0005]);
title('Difference of spectra at same points 100ms vs 200ms');
legend('Location', 'eastoutside');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
setSetting('plotName',mkNewDir(curSaveDir, strcat(altNames{k} ,'-pointSpectraDiff')));
savePlot(4);
    
figure(5);clf;
load('parameters/solax_reconstructed_spectrum.mat', 'solaxSpec', 'solaxLocalMaxWavelengths');
hold on 
h1 = plot(wavelengths, squeeze(mean(mean(spectVals(1,:,:,:), 2), 3)), 'DisplayName', 'Point Mean for 100ms', 'LineWidth', 5);
h2 = plot(wavelengths, squeeze(mean(mean(spectVals(2,:,:,:), 2), 3)), 'DisplayName', 'Point Mean for 200ms', 'LineWidth', 5);
h3 = plot(wavelengths, solaxSpec ./ 10^(6), 'DisplayName', 'Solax-iO illumination spectrum (div by 10^6)', 'LineWidth', 5);
for i = 1:numel(solaxLocalMaxWavelengths)
    xline(solaxLocalMaxWavelengths(i), '-', sprintf('Solax-iO Local Max at %d nm', solaxLocalMaxWavelengths(i)), 'LineWidth', 3);
end 
hold off 
ylim([0,0.0003]);
title('Spectra mean for selected spectra at same points 100ms vs 200ms');
legend([h1, h2, h3], 'Location', 'eastoutside');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
setSetting('plotName',mkNewDir(curSaveDir, strcat(altNames{k} ,'-pointSpectraMean')));
savePlot(5);