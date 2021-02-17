
%% Experiments
experiment = 'fusion_comparison';
%'effectOfFilter', 'effectOfAveraging', 'effectOfFusion',
%'integration_time_comparison', 'polarizing_effect_on_tissue'
% 'spatial_average_comparison', 'normalizationWithLR'

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
    case 'normalizationWithLR'
                
        experiment = 'normalizationWithLR';

        %% Normalize image to avoid spatial discrepancies using linera regression 
        dataDate = '20210127';
        configuration = 'singleLightClose';
        integrationTime = 1460;
        indirFolder = 'image fusion';
        % toFuseNames = {'20210127_122055_top_left_1', '20210127_122055_top_left_2', '20210127_122055_top_left_2b'};
        targetPosition = 'colorchartTopLeft';

        isNormalized = false; 
        hasWindow = true;
        plotSaveDir = fullfile(getSetting('savedir'), getSetting('saveFolder'));
        initialization;


        if isNormalized 
            setSetting('saveFolder', strcat(experiment, 'normalizedByPixel')); 
        else 
             setSetting('saveFolder', strcat(experiment, 'raw')); 
        end 

        if hasWindow 
            windowHalfSize = 300; 
            height = (windowHalfSize * 2 + 1);
            width = (windowHalfSize * 2 + 1);
            centerx = windowHalfSize+1;
            centery = centerx; 

        else 
            %Needs 12GB memory
            height = 1376;    
            width = 1024;
            centerx = floor(height/2) + 1;
            centery = floor(width/2) + 1;
        end
        k = 3;
        z = 401;
        xVals = zeros(k, height * width, z);
        yVals = zeros(k, z);

        for i = 1:k  
            imgName = strcat(targetPosition, num2str(i));

            filename = getFilename(getSetting('configuration'), imgName, getSetting('integrationTime'));
            [spectralData, ~, ~] = loadH5Data(filename, getSetting('experiment'));
            [m, n, ~] = size(spectralData);

            if hasWindow
                chartMask = zeros(m, n);
                chartMask( (floor(m/2) - windowHalfSize):(floor(m/2) + windowHalfSize), (floor(n/2) - windowHalfSize):(floor(n/2) + windowHalfSize) ) = 1;
                if isNormalized
                    croppedSpectra = readHSI(spectralData, chartMask);
                else
                    croppedSpectra = readHSI(spectralData, chartMask, 'raw');
                end
            else 
                croppedSpectra = spectralData; 
                clear 'spectralData';
            end 
            xVals(i, :, :) = reshape(croppedSpectra, [ height*width , 401]);
            yVals(i, :) = croppedSpectra(centerx, centery);
        end 

        %% Settings for plots 
        xVals4D = reshape(xVals, [k, height, width, z]);
        x1 = 10; 
        y1 = 10;
        x2 = 590;
        y2 = 590;
        pos = [x1 y1;x1 y2; x2 y1; x2 y2; centerx centery];
        startx = floor(m/2) - centerx + 1;
        starty = floor(n/2) - centery + 1;
        posXY(:,1) = pos(:,1) + startx;
        posXY(:,2) = pos(:,2) + starty;
        posXY = [posXY(:,2), posXY(:,1)]; % because data is reversed on the image
        colors = {'y', 'g', 'c', 'b', 'm'};
        curveNames = {'topleft', 'topright', 'bottomleft', 'bottomright', 'center'};

        %% Indexes for bandmax 
        [rref, inds] = max(reshape(croppedSpectra, [height * width, z]), [], 1);
        xy = zeros(z, 2);
        for i = 1:z
            [row, col] = ind2sub( [height, width], inds(i));
            xy(i, :) = [row; col];
        end 
        a = unique(inds');
        hc = [a, histc(inds(:),a)];
        freqInd = find(hc(:, 2) == max(hc(:,2)));
        [row, col] = ind2sub([height, width], a(freqInd));
        maxfreq = hc(freqInd, 2);
        fprintf('Pixel (%d, %d) has bandmax value %d times (most common) \n', col, row, maxfreq);

        fig = figure(3); clf;
        imshow(baseImage);
        hold on
        h = imshow(chartMask * 0.6) ;
        hold off
        set(h, 'AlphaData', chartMask);
        hold on 
        for i = 1:z
        scatter(xy(i,2) + starty, xy(i,1) + startx, 20, 'x', 'LineWidth',5,'MarkerEdgeColor','blue');
        end
        scatter(col + starty, row + startx, 20, 'x', 'LineWidth',5,'MarkerEdgeColor','red');
        hold off
        setSetting('plotName', fullfile( plotSaveDir, 'bandMaxPoints.png'));
        savePlot(fig);


        %% Linear Regression
        %lrReferencePos = [centerx, centery];
        lrReferencePos = [row, col];
        lrReferenceInd = sub2ind([height, width], row, col);
        yVals = xVals(:,lrReferenceInd, :);
        coeffs = zeros(height * width, z);
        for i = 1:height*width
            for j = 1:z
                coeffs(i, j) = xVals(:, i, j) \ yVals(:,j);
            end
        end

        v1 = coeffs;
        v2 = reshape(v1, [height, width, z]);


        fig = figure(1);clf;
        b = permute(v2(:,:,100), [2, 1]);
        imagesc(b);
        colorbar;
        title('At 480nm');
        setSetting('plotName', fullfile( plotSaveDir, 'lr_coeffs_480.png'));
        savePlot(fig);

        fig = figure(2);clf;
        b =  permute(v2(:,:,350), [2, 1]);
        imagesc(b);
        colorbar;
        title('At 730nm');
        setSetting('plotName', fullfile( plotSaveDir, 'lr_coeffs_730.png'));
        savePlot(fig);

        fig = figure(3);clf;
        v3 = permute(mean(v2, 3), [2,1]);
        imagesc(v3)
        colorbar;
        title('Average multiplier to adjust to center of screen level');
        setSetting('plotName', fullfile( plotSaveDir, 'lr_coeffs_average.png'));
        savePlot(fig);


        %% Show the mask that is being processed
        fig = figure(1); clf; 
        baseImage = rescale(spectralData(:,:,100));
        imagesc(baseImage);
        hold on
        h = imshow(chartMask * 0.6) ;
        hold off
        set(h, 'AlphaData', chartMask);
        hold on 
        for i = 1:numel(colors)
        scatter(posXY(i,1),posXY(i,2), 50, 'x', 'LineWidth',20,'MarkerEdgeColor',colors{i});
        end
        hold off 
        setSetting('plotName', fullfile( plotSaveDir, 'lr_mask.png'));
        savePlot(fig);


        %% Measured Spectra average for 5 masks 
        captureAvgWhite = zeros(numel(colors), z);
        for i = 1:numel(colors)
            captureAvgWhite(i,:) = mean(squeeze(xVals4D(:, pos(i,1), pos(i,2), :)), 1);
        end 
        plots(2, @plotColorChartSpectra, captureAvgWhite, curveNames, strcat('measured-raw', '_', '5masks'), ...
                    [0, 0.003], false);

        %% Measured Spectra 3 captures for 5 masks 
        captureAvgWhite = zeros(k, numel(colors), z);
        for i = 1:numel(colors)
            captureAvgWhite(:,i,:) = squeeze(xVals4D(:, pos(i,1), pos(i,2), :));
        end 
        captureAvgWhite = permute(captureAvgWhite, [2,3,1]);
        captureAvgWhitePlus = [captureAvgWhite; repmat(rref, [1, 1, 3])];
        curveNamesPlus = [curveNames, 'reference'];
        plots(4, @plotColorChartSpectra, captureAvgWhitePlus, curveNamesPlus, strcat('measured-raw', '_', '5masks3captures'), ...
                    [0, 0.003], false);

        %% After Normalization with max spectra 
        captureAvgWhite = zeros(numel(colors), z);
        for i = 1:numel(colors)
            captureAvgWhite(i,:) = squeeze(xVals4D(1, pos(i,1), pos(i,2), :))' ./ rref;
        end 
        plots(5, @plotColorChartSpectra, captureAvgWhite, curveNames, strcat('measured', '_', '5masksNormalized'), ...
                    [0, 1.4], false);

        %% After Multiplication with lr coeff and normalization with bandmax 
        captureAvgWhite = zeros(numel(colors), z);
        for i = 1:numel(colors)
            captureAvgWhite(i,:) = squeeze(xVals4D(1, pos(i,1), pos(i,2), :))' .* squeeze(v2(pos(i,1), pos(i,2), :))' ./ rref;
        end 
        plots(6, @plotColorChartSpectra, captureAvgWhite, curveNames, strcat('measured', '_', '5masksLRNormalized'), ...
                    [0, 1.4], false);



    case 'effectOfFilter' %'polarizing_filter_use_comparison'
        
        %% Compare effect of filter usage
        dataDate = '20210111';
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
        
    case 'effectOfAveraging' %'capture_average_comparison'
        
        %% Compare image averaging
        dataDate = '20210127';
        configuration = 'singleLightClose';
        integrationTime = 1460;
        indirFolder = 'image fusion';
%         toFuseNames = {'20210127_122055_top_left_1', '20210127_122055_top_left_2', '20210127_122055_top_left_2b'};
        targetPosition = 'colorchartTopLeft';
        initialization;
        
        setSetting('saveFolder', strcat(experiment, '\singlePixelMask')); %'\squarePatchMask'
        
        n = 3;
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
        
        plots(2, @plotColorChartSpectra, spectVals, curveNames, strcat('measured', '_', 'allResults'), ...
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
        
    case 'fusion_comparison' 
        % Compare image fusion
        dataDate = '20210127';
        configuration = 'singleLightClose';
        integrationTime = 1460;
        indirFolder = 'image fusion';
        targetPosition = 'colorchartTopLeft'; %'colorchartTopLeft' 'colorchartBottomRight''colorchartBottomLeft'
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
        
        %% For light skin only 
        plots(1, @plotColorChartSpectra, squeeze(multipleSpectra(2,:,:))', {'Fused1', 'Fused2', 'Fused3', 'Fused4', 'Fused5', 'Fused6'}, 'measured_fusedLightSkin', [], false);
        
    case 'integration_time_comparison' %'integration_time_comparison'
        
        %% Compare different exposure times
        dataDate = '20201209';
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
        
    case 'polarizing_effect_on_tissue' %'polarizing_effect_on_tissue'
        
        %% Compare effect of polarizing on skin
        dataDate = '20210127';
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
        
    case 'spatial_average_comparison' %'spatial_average_comparison'
        
        %% Compare effect of spatial averaging
        dataDate = '20201220';
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