
%% Experiments
experiment = 'effectOfFusion';
%'effectOfFilter', 'effectOfAveraging', 'effectOfFusion',
%'integration_time_comparison', 'polarizing_effect_on_tissue'
% 'spatial_average_comparison'

%% Read white
readBlackWhite = false;
if readBlackWhite
    for i = 1:2
        setSetting('saveFolder', fullfile(experiment, configurations{i}));
        readWhite(dataDate, integrationTime, true, false, experiment, configuration, 'white', [], hasFilter{i});
    end
    
    for i = [500, 800, 1460] %1800,
        integrationTime = i;
        setSetting('integrationTime', integrationTime);
        readWhite(dataDate, integrationTime, true, false, experiment, configuration, 'white');
    end
end

switch experiment
    case 'effectOfFilter'
        
        %% Compare effect of filter usage
        dataDate = '20210111';
        experiment = 'polarizing_filter_use_comparison';
        integrationTime = 1460;
        initialization;
        
        % Available configurations 'no_filter', 'filter'
        configurations = {'no_filter', 'filter'};
        hasFilter = {false, true};
        
        n = 8;
        tables = cell(n, 1);
        measuredSpectra = cell(n, 1);
        adjustedSpectra = cell(n, 1);
        allowRoiSelection = true;
        ind = 0;
        for i = 1:2
            setSetting('saveFolder', fullfile(experiment, configurations{i}));
            
            datafiles = dir(fullfile(getSetting('datadir'), '*.h5'));
            whiteFilename = datafiles(contains({datafiles.name}, 'white')).name;
            getRepresentativePoints(whiteFilename);
            
            positions = {'left_up', 'left_down', 'right_down', 'right_up'};
            
            for j = 1:4
                close all;
                position = positions{j};
                getSetting('targetPosition', position);
                positionFilename = datafiles(contains({datafiles.name}, position)).name;
                setSetting('saveFolder', strcat(experiment, '\', position, '_', configurations{i}));
                ind = ind + 1;
                [tables{ind}, measuredSpectra{ind}, adjustedSpectra{ind}] = evaluateColorchart(positionFilename, allowRoiSelection);
            end
        end
        
    case 'effectOfAveraging'
        
        %% Compare image averaging
        dataDate = '20210127';
        experiment = 'capture_average_comparison';
        configuration = 'singleLightClose';
        integrationTime = 1460;
        indirFolder = 'image fusion';
        toFuseNames = {'20210127_122055_top_left_1.h5', '20210127_122055_top_left_2.h5', '20210127_122055_top_left_2b.h5'};
        targetPosition = 'colorchartTopLeft';
        initialization;
        
        n = length(toFuseNames);
        
        setSetting('saveFolder', strcat(experiment, '\singlePixelMask')); %'\squarePatchMask'
        %         curSaveDir = mkNewDir(getSetting('savedir'), getSetting('saveFolder'));
        
        xLen = 6;
        yLen = 5;
        z = 401;
        spectVals = zeros(n, xLen, yLen, z);
        
        for k = 1:n
            imgName = strcat(targetPosition, num2str(k));
            setSetting('saveFolder', fullfile(experiment, targetPosition, num2str(k)));
            [spectVals(k, :, :, :), curveNames] = getRepresentativePoints(imgName, 50:200:1088, 50:200:982);
        end
        
        limits = [0, 0.008];
        spectValsAverage = squeeze(mean(spectVals, 1));
        plots(1, @plotColorChartSpectra, spectValsAverage, curveNames, strcat('measured', '_', imgName, '3captureAverage'), ...
            limits, false);
        
        plots(2, @plotColorChartSpectra, spectVals, curveNames, strcat(curCase, '_', 'allResults'), ...
            limits, true);
        
        % Compare curves
        close all;
        tables = cell(n+1, 1);
        measuredSpectra = cell(n+1, 1);
        adjustedSpectra = cell(n+1, 1);
        allowRoiSelection = true;
        avgMeasured = zeros(24, 36);
        
        for k = 1:n
            imgName = strcat(targetPosition, num2str(k));
            filename = getFilename(configuration, imgName);
            setSetting('saveFolder', fullfile(experiment, targetPosition, num2str(k)));
            [tables{k}, measuredSpectra{k}, adjustedSpectra{k}] = evaluateColorchart(filename, allowRoiSelection);
            avgMeasured = avgMeasured + measuredSpectra{k};
        end
        measuredSpectra{n+1} = avgMeasured ./ n;
        expected = getExpectedValues();
        tables{n+1} = compareSpectra(expected, measuredSpectra{n+1}, [tables{k}.Patch]);
        
        for i = 1:(n + 1)
            tables{i, 1}{25, 1} = {'Average'};
            tables{i, 1}{26, 1} = {'StandardDeviation'};
            for j = 2:7
                tables{i, 1}{25, j} = mean(tables{i, 1}{:, j});
                tables{i, 1}{26, j} = std(tables{i, 1}{:, j});
            end
        end
        
        for j = 1:24
            if tables{1, 1}.AdjGoF(j) < 0.8
                fprintf('**Low result for patch %s\n', tables{1, 1}.Patch{j});
            end
        end
        
        savedir = fullfile(getSetting('savedir'), getSetting('saveFolder'), 'resultTable.mat');
        save(savedir, 'tables', 'measuredSpectra', 'adjustedSpectra');
        fprintf('Saved values in %s,\n', savedir);
        
        for i = 1:(n + 1)
            writetable(tables{i, 1}, fullfile(getSetting('savedir'), getSetting('saveFolder'), 'resultTable.xlsx'), 'Sheet', i);
        end
        
        % For measured
        close all;
        multipleSpectra = zeros(size(measuredSpectra{1}, 1), size(measuredSpectra{1}, 2), n);
        for i = 1:n
            multipleSpectra(:, :, i) = measuredSpectra{i};
        end
        [~, patchNames] = getExpectedValues();
        plots(1, @plotColorChartSpectra, multipleSpectra, patchNames, 'measured_imageAverage', [], true);
        
        % For adjusted
        close all;
        multipleSpectra = zeros(size(adjustedSpectra{1}, 1), size(adjustedSpectra{1}, 2), n);
        for i = 1:n
            multipleSpectra(:, :, i) = adjustedSpectra{i};
        end
        plots(1, @plotColorChartSpectra, multipleSpectra, patchNames, 'measured-adjusted_imageAverage', [], true);
        
    case 'effectOfFusion'
        % Compare image fusion
        dataDate = '20210127';
        experiment = 'fusion_comparison';
        configuration = 'singleLightClose';
        integrationTime = 1460;
        indirFolder = 'image fusion';
        targetPosition = 'colorchartTopRight'; %'colorchartTopLeft' 'colorchartBottomRight''colorchartBottomLeft'
        initialization;
        n = 8;
        
        % Select various points on the image and visualize spectra
        xLen = 6;
        yLen = 5;
        z = 401;
        spectVals = zeros(n, xLen, yLen, z);
        for k = 1:n
            if k ~= 3
                imgName = strcat(targetPosition, num2str(k));
                setSetting('saveFolder', fullfile(experiment, targetPosition, num2str(k)));
                [spectVals(k, :, :, :), ~] = getRepresentativePoints(imgName, 50:200:1088, 50:200:982);
            end
        end
        
        %% Compare colorchart curves per single image
        close all;
        tables = cell(n+1, 1);
        measuredSpectra = cell(n+1, 1);
        adjustedSpectra = cell(n+1, 1);
        allowRoiSelection = true;
        
        ind = 0;
        for k = 1:9
            if k ~= 3
                ind = ind + 1;
                imgName = strcat(targetPosition, num2str(k));
                [filename, integrationTime] = getFilename(configuration, imgName);
                setSetting('saveFolder', fullfile(experiment, targetPosition, num2str(k)));
                [tables{ind}, measuredSpectra{ind}, adjustedSpectra{ind}] = evaluateColorchart(filename, allowRoiSelection);
            end
        end
        
        %% Fusion
        fusedSpectra = cell(6, 1);
        
        % Single image at 1460ms
        avgMeasured1 = measuredSpectra{1};
        fusedSpectra{1} = avgMeasured1;
        
        % Fuse 2 images at 1460ms
        avgMeasured2 = (measuredSpectra{1} + measuredSpectra{2}) / 2;
        fusedSpectra{2} = avgMeasured2;
        
        % Single image at 1800ms
        avgMeasured3 = measuredSpectra{3};
        fusedSpectra{3} = avgMeasured3;
        
        % Fuse 2 images at 1460 (380-540nm) and 1800 (541-780nm)
        [~, ~, ia] = intersect(getWavelengths(36), [380:540]);
        
        avgMeasured4 = zeros(24, 36);
        avgMeasured4(:, 1:length(ia)) = measuredSpectra{4};
        avgMeasured4(:, length(ia)+1:end) = measuredSpectra{5};
        fusedSpectra{4} = avgMeasured4;
        
        % Fusion image at 1460 and 1800
        avgMeasured5 = zeros(24, 36);
        avgMeasured5(:, 1:length(ia)) = measuredSpectra{1}(:, 1:length(ia));
        avgMeasured5(:, length(ia)+1:end) = (measuredSpectra{1}(:, length(ia) + 1:end) + measuredSpectra{5}) / 2;
        fusedSpectra{5} = avgMeasured5;
        
        % Fusion image at 2x500 ms (380-780nm) and 1x800ms (541-780nm)
        avgMeasured6 = (measuredSpectra{6} + measuredSpectra{7}) / 2;
        avgMeasured6(:, length(ia)+1:end) = (avgMeasured6(:, length(ia)+1:end) + measuredSpectra{5}) / 2;
        fusedSpectra{6} = avgMeasured6;
        
        % Compare colorchart Curves per fused image
        n = length(fusedSpectra);
        fusedSpectraTables = cell(n, 1);
        [expected, patchNames] = getExpectedValues();
        
        for i = 1:n
            setSetting('saveFolder', fullfile(experiment, targetPosition, strcat('fusedSpectra', num2str(i))));
            fusedSpectraTables{i} = compareSpectra(expected, fusedSpectra{i}, [tables{1}.Patch]);
        end
        
        for i = 1:n
            fusedSpectraTables{i, 1}{25, 1} = {'Average'};
            fusedSpectraTables{i, 1}{26, 1} = {'StandardDeviation'};
            for j = 2:7
                fusedSpectraTables{i, 1}{25, j} = mean(fusedSpectraTables{i, 1}{:, j});
                fusedSpectraTables{i, 1}{26, j} = std(fusedSpectraTables{i, 1}{:, j});
            end
        end
        
        %% Show failed results
        for j = 1:24
            if fusedSpectraTables{1, 1}.AdjGoF(j) < 0.8
                fprintf('**Low result for patch %s\n', fusedSpectraTables{1, 1}.Patch{j});
            end
        end
        
        setSetting('saveFolder', fullfile(experiment, targetPosition));
        savedir = fullfile(getSetting('savedir'), getSetting('saveFolder'), 'resultTable.mat');
        save(savedir, 'tables', 'measuredSpectra', 'adjustedSpectra', 'fusedSpectra', 'fusedSpectraTables');
        fprintf('Saved values in %s,\n', savedir);
        
        for i = 1:n
            writetable(fusedSpectraTables{i, 1}, fullfile(getSetting('savedir'), getSetting('saveFolder'), 'resultTable.xlsx'), 'Sheet', i);
        end
        
        %% Comparison for all fused spectra
        close all;
        multipleSpectra = zeros(size(fusedSpectra{1}, 1), size(fusedSpectra{1}, 2), n);
        for i = 1:n
            multipleSpectra(:, :, i) = fusedSpectra{i};
        end
        plots(1, @plotColorChartSpectra, multipleSpectra, patchNames, 'measured_fused', [], false);
        
    case 'integration_time_comparison'
        
        %% Compare different exposure times
        dataDate = '20201209';
        experiment = 'integration_time_comparison';
        configuration = 'singleLightClose';
        integrationTime = 1460;
        initialization;
        
        setSetting('saveFolder', experiment);
        curSaveDir = mkNewDir(getSetting('savedir'), getSetting('saveFolder'));
        exp100name = '20201209_173035-manualexposure100.h5';
        exp200name = '20201209_171535-manualexposure200.h5';
        names = {exp100name, exp200name};
        altNames = {'int100', 'int200'};
        
        xLen = 7;
        yLen = 7;
        z = 401;
        spectVals = zeros(2, xLen, yLen, z);
        for k = 1:2
            imgName = names{k};
            setSetting('saveFolder', fullfile(experiment, altNames{k}));
            [spectVals(k, :, :, :), curveNames] = getRepresentativePoints(imgName, 50:200:size(spectralData, 1), 50:200:size(spectralData, 2));
        end
        
        %%Difference of spectra at same points 100ms vs 200ms
        spectra1 = reshape(squeeze(spectVals(1, i, j, :)), [xLen * yLen, z]);
        spectra2 = reshape(squeeze(spectVals(2, i, j, :)), [xLen * yLen, z]);
        diffSpect = reshape(spectra1-spectra2, [xLen, yLen, z]);
        limits = [-0.0005, 0.0005];
        plots(1, @plotColorChartSpectra, diffSpect, curveNames, strcat('difference', '_', 'pointSpectraDiff'), ...
            limits, false);
        
        
        figure(5);
        clf;
        load('parameters/solax_reconstructed_spectrum.mat', 'solaxSpec', 'solaxLocalMaxWavelengths');
        hold on
        h1 = plot(wavelengths, squeeze(mean(mean(spectVals(1, :, :, :), 2), 3)), 'DisplayName', 'Point Mean for 100ms', 'LineWidth', 5);
        h2 = plot(wavelengths, squeeze(mean(mean(spectVals(2, :, :, :), 2), 3)), 'DisplayName', 'Point Mean for 200ms', 'LineWidth', 5);
        h3 = plot(wavelengths, solaxSpec./10^(6), 'DisplayName', 'Solax-iO illumination spectrum (div by 10^6)', 'LineWidth', 5);
        for i = 1:numel(solaxLocalMaxWavelengths)
            xline(solaxLocalMaxWavelengths(i), '-', sprintf('Solax-iO Local Max at %d nm', solaxLocalMaxWavelengths(i)), 'LineWidth', 3);
        end
        hold off
        ylim([0, 0.0003]);
        title('Spectra mean for selected spectra at same points 100ms vs 200ms');
        legend([h1, h2, h3], 'Location', 'eastoutside');
        set(gcf, 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);
        setSetting('saveFolder', experiment);
        setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('saveFolder'), strcat(altNames{k}, '-pointSpectraMean')));
        savePlot(5);
        
    case 'polarizing_effect_on_tissue'
        
        %% Compare effect of polarizing on skin
        dataDate = '20210127';
        experiment = 'polarizing_effect_on_tissue';
        configuration = 'singleLightClose';
        indirFolder = 'skin test';
        initialization;
        
        imageFilenames = {'20210127_174840_dry_nofilter_800.h5', ...
            '20210127_174840_wet_filter_1460.h5', ...
            '20210127_174840_dry_filter_1460.h5', ...
            '20210127_174840_wet_nofilter_300.h5', ...
            '20210127_174840_dry_nofilter_1460.h5', ...
            '20210127_174840_dry_nofilter_300.h5'};
        integrationTimes = [800, 1460, 1460, 300, 1460, 300];
        dryInd = [3, 6];
        wetInd = [2, 4];
        
        setSetting('saveFolder', experiment);
        
        xPoints = 300:150:800;
        yPoints = 300:150:800;
        allowRoiSelection = true;
        limits = {[0, 0.0035], [0, 0.015]};
        actualSpectralVals = cell(4, 1);
        
        %% For dry skin
        baseName = 'dryHand';
        
        for k = 1:length(dryInd)
            ind = dryInd(k);
            hasFilter = ~contains(imageFilenames{ind}, 'nofilter');
            if hasFilter
                imgName = strcat(baseName, 'Filter');
            else
                imgName = strcat(baseName, 'NoFilter');
            end
            integrationTime = integrationTimes(ind);
            setSetting('integrationTime', integrationTime);
            
            setSetting('saveFolder', fullfile(experiment, baseName, altNames{k}));
            getRepresentativePoints(imgName, xPoints, yPoints);
            
            filename = getFilename(getSetting('configuration'), imgName, integrationTime);
            [spectralData, ~, ~] = loadH5Data(filename, experiment);
            [colorMasks, chartMask] = getColorchartMasks(squeeze(spectralData(:, :, 100)), allowRoiSelection, experiment); %(imageXYZ, allowRoiSelection, experiment)
            actualSpectralVals{k} = readHSI(spectralData, {chartMask, colorMasks}, [], hasFilter);
            lineNames = arrayfun(@num2str, [1:30], 'UniformOutput', false);
            plots(1, @plotColorChartSpectra, actualSpectralVals{k}, lineNames, strcat('measured_', baseName, altNames{k}), limits{1});
        end
        
        %% For wet skin
        close all;
        baseName = 'wetHand';
        
        for k = 1:length(wetInd)
            ind = wetInd(k);
            hasFilter = ~contains(imageFilenames{ind}, 'nofilter');
            if hasFilter
                imgName = strcat(baseName, 'Filter');
            else
                imgName = strcat(baseName, 'NoFilter');
            end
            integrationTime = integrationTimes(ind);
            setSetting('integrationTime', integrationTime);
            
            setSetting('saveFolder', fullfile(experiment, baseName, altNames{k}));
            getRepresentativePoints(imgName, xPoints, yPoints);
            
            filename = getFilename(getSetting('configuration'), imgName, integrationTime);
            [spectralData, ~, ~] = loadH5Data(filename, experiment);
            [colorMasks, chartMask] = getColorchartMasks(squeeze(spectralData(:, :, 100)), allowRoiSelection, experiment); %(imageXYZ, allowRoiSelection, experiment)
            actualSpectralVals{k} = readHSI(spectralData, {chartMask, colorMasks}, [], hasFilter);
            lineNames = arrayfun(@num2str, [1:30], 'UniformOutput', false);
            plots(1, @plotColorChartSpectra, actualSpectralVals{k}, lineNames, strcat('measured_', baseName, altNames{k}), limits{2});
        end
        
        setSetting('saveFolder', experiment);
        savedir = fullfile(getSetting('savedir'), getSetting('saveFolder'), 'normalizedSpectra.mat');
        save(savedir, 'actualSpectralVals');
        fprintf('Save actual spectral values in %s\n\n', savedir);
        
        %% Show images with error
        
        setSetting('saveImages', false);
        close all;
        %800
        setSetting('integrationTime', 800);
        imgName = 'dryHandNoFilter';
        getRepresentativePoints(imgName, xPoints, yPoints, [0, 0.01]);
        
        %1460
        setSetting('integrationTime', 1460);
        imgName = 'dryHandNoFilter';
        getRepresentativePoints(imgName, xPoints, yPoints, [0, 0.01]);
        
    case 'spatial_average_comparison'
        
        %% Compare effect of spatial averaging
        dataDate = '20201220';
        experiment = 'spatial_average_comparison';
        configuration = 'doubleLightClose';
        initialization;
        
        setSetting('saveFolder', experiment);
        
        curSaveDir = mkNewDir(getSetting('savedir'), getSetting('saveFolder'));
        
        names = {'doubleLightClose', 'singleLightClose', 'singleLightFar'};
        spectVals = zeros(numel(names), 6, 5, 401);
        spect2Vals = zeros(numel(names), 6, 5, 401);
        spect3Vals = zeros(numel(names), 6, 5, 401);
        
        xPoints = 50:200:1088;
        yPoints = 50:200:982;
        
        for k = 1:numel(names)
            filename = getFilename(names{k}, 'colorchart');
            tempDataDate = strsplit(filename, '_');
            indir = strcat(fullfile(originDir, '2_saitamaHSI\'), tempDataDate{1}, '_test\h5\');
            setSetting('datadir', indir);
            [spectVals(k, :, :, :), curveNames] = getRepresentativePoints(imgName, xPoints, yPoints, [], 1);
            spect2Vals(k, :, :, :) = getRepresentativePoints(imgName, xPoints, yPoints, [], 2);
            spect3Vals(k, :, :, :) = getRepresentativePoints(imgName, xPoints, yPoints, [], 3);
        end
        
        load('parameters/solax_reconstructed_spectrum.mat', 'solaxSpec', 'solaxLocalMaxWavelengths');
        nn = numel(xPoints);
        mm = numel(yPoints);
        allCurves = zeros(3, 3, 401);
        allCurves(3, 1, :) = removePoints(spectVals(3, :, :, :), nn, mm, [1:10, 18, 19, 23, 24]);
        allCurves(3, 2, :) = removePoints(spect2Vals(3, :, :, :), nn, mm, [1:10, 18, 19, 23, 24]);
        allCurves(3, 3, :) = removePoints(spect3Vals(3, :, :, :), nn, mm, [1:10, 18, 19, 23, 24]);
        allCurves(2, 1, :) = removePoints(spectVals(2, :, :, :), nn, mm, [13, 18]);
        allCurves(2, 2, :) = removePoints(spect2Vals(2, :, :, :), nn, mm, [13, 18]);
        allCurves(2, 3, :) = removePoints(spect3Vals(2, :, :, :), nn, mm, [13, 18]);
        allCurves(1, 1, :) = removePoints(spectVals(1, :, :, :), nn, mm, [13, 18]);
        allCurves(1, 2, :) = removePoints(spect2Vals(1, :, :, :), nn, mm, [13, 18]);
        allCurves(1, 3, :) = removePoints(spect3Vals(1, :, :, :), nn, mm, [13, 18]);
        
        figure(5);
        clf;
        hold on
        k = 0;
        h = zeros(9, 1);
        for i = 1:3
            for j = 1:3
                k = k + 1;
                if j == 1
                    mark = '-';
                elseif j == 2
                    mark = '--';
                else
                    mark = '.';
                end
                dispName = sprintf('%s + %dx%d', names{i}, j, j);
                h(k) = plot(wavelengths, squeeze(allCurves(i, j, :)), mark, 'DisplayName', dispName, 'LineWidth', 3);
            end
        end
        h(10) = plot(wavelengths, solaxSpec./10^(5), 'DisplayName', 'Solax-iO illumination spectrum (div by 10^5)', 'LineWidth', 5);
        for i = 1:numel(solaxLocalMaxWavelengths)
            xline(solaxLocalMaxWavelengths(i), '-', sprintf('Solax-iO Local Max at %d nm', solaxLocalMaxWavelengths(i)), 'LineWidth', 3);
        end
        hold off
        ylim([0, 0.005]);
        title('Spectra mean for selected spectra with spatial averaging');
        legend(h, 'Location', 'eastoutside');
        set(gcf, 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);
        setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('saveForlder'), strcat(experiment, '-pointSpectraMean')));
        savePlot(5);
        
    otherwise
        error('Unsupported experiment')
end

%% Function for 'spatial_average_comparison'
function [res] = removePoints(vals, nn, mm, toRemove)
valsTemp = reshape(squeeze(vals), [nn * mm, 401]);
valsTemp2 = zeros(nn*mm-numel(toRemove), 401);
j = 0;
for i = 1:size(valsTemp, 1)
    if isempty(find(toRemove == i, 1))
        j = j + 1;
        valsTemp2(j, :) = valsTemp(i, :);
    end
end
res = squeeze(mean(valsTemp2, 1));
end