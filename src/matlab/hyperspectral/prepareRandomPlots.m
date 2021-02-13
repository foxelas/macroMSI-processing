%% Prepare graphs 
close all;

%% What to run 
runsWhiteRef = false; 
runsWhitePatch = true; 


%% Settings 
baseName = 'colorchartBottomLeft'; %'colorchartTopLeft' 'colorchartBottomRight''colorchartBottomLeft' 'colorchartTopRight'
experiment = 'random_graphs';
dataDate = '20210127';
initialization; 

%% Get white sheet values
if runsWhiteRef 
    maskN = 5;
    masks = cell(maskN,1);
    for i = [ 500, 800, 1460, 1800]
        integrationTime = i;
        setSetting('integrationTime', integrationTime);
        loadDir =  fullfile( matdir, configuration, strcat(num2str(integrationTime), '_fullreflectance_ByPixel.mat'));
        load(loadDir);
        [m,n,l] = size(fullReflectanceByPixel);
        spectralData = fullReflectanceByPixel;
        baseImage = rescale(spectralData(:,:,300));

        if (i == 500)
            mask = zeros(m,n);
            mask(50:60, 50:60) = 1;
            masks{1} = mask; 
            mask = zeros(m,n);
            mask(50:60, 800:810) = 1;
            masks{2} = mask;
            mask = zeros(m,n);
            mask(1320:1330, 50:60) = 1;
            masks{3} = mask;
            mask = zeros(m,n);
            mask(1320:1330, 800:810) = 1;
            masks{4} = mask;
            mask = zeros(m,n);
            mask(700:710, 500:510) = 1;
            masks{5} = mask;

            fiveMasks = (masks{1} + masks{2} + masks{3} + masks{4} + masks{5}) / 5; 

            position =  [ 55 55 ; 55 1325; 805 55; 805 1325; 505 705];
            roiNames = {'mask1', 'mask2', 'mask3', 'mask4', 'mask5'};
            annotatedBaseImage = insertText(baseImage,position,roiNames,'AnchorPoint','LeftBottom', 'FontSize', 20);

            fig2 = figure(2);
            s = imshow(annotatedBaseImage);
            alpha(s,.5);
            hold on;
            h = imshow(fiveMasks);
            hold off;
            set(h, 'AlphaData', fiveMasks);
            setSetting('saveFolder', strcat(experiment,'\'));

            setSetting('plotName',mkNewDir(getSetting('savedir'), getSetting('saveFolder'), 'masksForReferenceSheet'));
            savePlot(fig2);
        end 

        averagePerMask = zeros(maskN, 401, 1);
        for j = 1:maskN 
            croppedSpectra = readHSI(spectralData, masks{j}, 'raw');
            croppedSpectra = reshape(croppedSpectra, [size(croppedSpectra, 1) * size(croppedSpectra,2), size(croppedSpectra, 3)]);
            averagePerMask(j, :) = mean(croppedSpectra);
        end
        plots(1, @plotColorChartSpectra, averagePerMask, roiNames, strcat('measured-raw_referenceRoiIntegration', num2str(integrationTime)), [0, 0.005]);

    end 
end

%% Prepare graphs for white patch before normalization 

if runsWhitePatch
    close all; 
    n = length(8);
    tables = cell(n+1, 1);
    measuredSpectra = cell(n+1,1);
    adjustedSpectra = cell(n+1,1);
    allowRoiSelection = true;

    curveNames = {'at 1460ms (a)', 'at 1460ms (b)', 'at 1800ms', 'at 1460ms (low)', ...
        'at 1800ms (high)', 'at 500ms (a)', 'at 500ms (b)', 'at 800ms (high)' };
    
    whitePatchRaw = zeros(8, 401);
    whitePatchNormalized = zeros(8, 401);
    lightSkinPatchNormalized = zeros(8, 401);

    ind = 0;  
    for k = 1:9
        if k ~= 3 
            ind = ind + 1;
            imgName = strcat(baseName, num2str(k));
            [filename, integrationTime]= getFilename(configuration, imgName);
            indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'),'image fusion', 'h5');
            setSetting('datadir', indir);    
            indexWhite = 8;
            indexLightSkin = 2; 
            patchIndexes = [indexLightSkin, indexWhite];
            
            [~, measuredSpectraRaw, ~] = evaluateColorchart(filename, experiment, allowRoiSelection, patchIndexes, 'raw');
            whitePatchRaw(ind, :) = pads(squeeze(measuredSpectraRaw(:, 2)));
            
            [~, measuredSpectra, ~] = evaluateColorchart(filename, experiment, allowRoiSelection, patchIndexes);
            lightSkinPatchNormalized(ind, 1:length(croppedSpectra)) = pads(squeeze(measuredSpectra(:,1)));
            whitePatchNormalized(ind, 1:length(croppedSpectra)) = pads(squeeze(measuredSpectra(:,2)));
        end 
    end
    
    setSetting('saveFolder', strcat(experiment,'\', baseName));
    plots(2, @plotColorChartSpectra, whitePatchRaw, curveNames, 'measured-raw_whitePatch', [0, 0.003]);
    plots(3, @plotColorChartSpectra, whitePatchNormalized, curveNames, 'measured-raw_whitePatchNormalized', [0,2]);
    plots(4, @plotColorChartSpectra, lightSkinPatchNormalized, curveNames, 'measured-raw_whitePatchNormalized', [0,2]);

end 