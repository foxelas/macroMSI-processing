%% Experiments 
experiment = 'effectOfFilter'; 

% Read white 
readBlackWhite = false;
if readBlackWhite
    for i=1:2
        setSetting('saveFolder', fullfile(experiment, configurations{i}));
            readWhite(dataDate, integrationTime, true, false, experiment, configuration, 'white', [], hasFilter{i} );
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
        measuredSpectra = cell(n,1);
        adjustedSpectra = cell(n,1);
        allowRoiSelection = true;
        ind = 0;
        for i = 1:2 
            setSetting('saveFolder', fullfile(experiment, configurations{i}));

            datafiles = dir(fullfile(getSetting('datadir'), '*.h5'));
            whiteFilename = datafiles( contains({datafiles.name}, 'white')).name;
            getRepresentativePoints(whiteFilename);

            positions = {'left_up', 'left_down', 'right_down', 'right_up'};

            for j = 1:4
                close all; 
                position = positions{j};

                positionFilename = datafiles( contains({datafiles.name}, position)).name;
                setSetting('saveFolder', strcat(experiment, '\', position, '_', configurations{i}));
                ind = ind + 1;
                [tables{ind}, measuredSpectra{ind}, adjustedSpectra{ind}] = evaluateColorchart(positionFilename, experiment, allowRoiSelection);    
            end 
        end 


    case ''
    %% Compare image averaging 
        dataDate = '20210127';
        experiment = 'capture_average_comparison';
        configuration = 'singleLightClose';
        integrationTime = 1460;
        toFuseNames = {'20210127_122055_top_left_1.h5', '20210127_122055_top_left_2.h5', '20210127_122055_top_left_2b.h5'};
        baseName = 'colorchartTopLeft';
        initialization;

        n = length(toFuseNames);

        setSetting('saveFolder', strcat(experiment, '\singlePixelMask')); %'\squarePatchMask'
%         curSaveDir = mkNewDir(getSetting('savedir'), getSetting('saveFolder'));

        spectVals =  zeros(n, 6, 5, 401);
        for k = 1:n
            imgName = strcat(baseName, num2str(k));
            filename = getFilename(configuration, imgName);
            
            [spectVals(k, :, :, :), curveNames ]= getRepresentativePoints(filename);
        end
            
        spectValsAverage = squeeze(mean(spectVals, 1));
        limits = [0,0.008];
        showMultiplePointSpectra(4, 1:6, 1:5, 0, getWavelengths(401), spectValsAverage, limits, strcat(imgName, '3captureAverage'));

        plots(5, @plotColorChartSpectra, spectVals, curveNames, strcat(curCase, '_', 'allResults'), ...
        limits, true);
      
%%%% Continue here 
        %% Compare curves 
        close all; 
        tables = cell(n+1, 1);
        measuredSpectra = cell(n+1,1);
        adjustedSpectra = cell(n+1,1);
        allowRoiSelection = true;
        avgMeasured = zeros(24, 36);

        for k = 1:n

            imgName = strcat(baseName, num2str(k));
            filename = getFilename(configuration, imgName);
            indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'),'image fusion', 'h5');
            setSetting('datadir', indir);

            [tables{k}, measuredSpectra{k}, adjustedSpectra{k}] = evaluateColorchart(filename, experiment, allowRoiSelection);
            avgMeasured = avgMeasured + measuredSpectra{k};
        end
        measuredSpectra{n+1} = avgMeasured ./ n; 
        expected = getExpectedValues();
        tables{n+1} = compareSpectra(expected, measuredSpectra{n+1}, [tables{k}.Patch]);

        for i = 1:(n+1)
        tables{i,1}{25,1} = {'Average'};
        tables{i,1}{26,1} = {'StandardDeviation'};
            for j = 2:7
            tables{i,1}{25,j} = mean(tables{i,1}{:,j});
            tables{i,1}{26,j} = std(tables{i,1}{:,j});
            end
        end 

        for j = 1:24
            if tables{1,1}.AdjGoF(j) < 0.8
                fprintf('**Low result for patch %s\n', tables{1,1}.Patch{j});
            end
        end

        savedir = fullfile(getSetting('savedir'), getSetting('saveFolder'), 'resultTable.mat');
        save(savedir, 'tables', 'measuredSpectra', 'adjustedSpectra');
        fprintf('Saved values in %s,\n', savedir);

        for i = 1:(n+1)
            writetable(tables{i,1}, fullfile(getSetting('savedir'), getSetting('saveFolder'), 'resultTable.xlsx'),'Sheet',i);
        end

        %% For measured 
        close all; 
        multipleSpectra = zeros(size(measuredSpectra{1},1), size(measuredSpectra{1},2), n);
        for i = 1:n
            multipleSpectra(:,:,i) = measuredSpectra{i};
        end 
        [~, patchNames] = getExpectedValues();
        plots(1, @plotColorChartSpectra, multipleSpectra, patchNames, 'measured_imageAverage');

        %% For adjusted 
        close all; 
        multipleSpectra = zeros(size(adjustedSpectra{1},1), size(adjustedSpectra{1},2), n);
        for i = 1:n
            multipleSpectra(:,:,i) = adjustedSpectra{i};
        end 
        plots(1, @plotColorChartSpectra, multipleSpectra, patchNames,'measured-adjusted_imageAverage');

end 