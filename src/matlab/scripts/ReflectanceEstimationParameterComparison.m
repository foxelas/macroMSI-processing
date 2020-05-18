if (false)

    %% Run comparison
    sigma1range = 10.^(-5:0.2:0);
    sigma2range = 10.^(-5:0.2:0);
    olympusSigmaRange = kron(10.^(-4:-3), [1, 5]);
    sigma1n = length(sigma1range);
    sigma2n = length(sigma2range);
    aRange = 0:0.2:1;
    an = length(aRange);
    gammaRange = 0.5:0.5:3;
    powerRange = -1:1;
    gamman = length(gammaRange);
    rhoRange = [0.2, 0.3, 0.4, 0.5, 0.6]; %[0.6, 0.65, 0.7, 0.8, 0.9, 0.97]
    windowRange = [3, 5, 9];

    setSetting('pixelValueSelectionMethod', 'extended');

    idx = 0;
    for rho = [0.6, 0.65, 0.7, 0.8, 0.9, 0.97]
        for windowDim = [3, 5, 9]
            for sigma = olympusSigmaRange
                idx = idx + 1;
            end
        end
    end
    for i = sigma1range
        idx = idx + 1;
    end
    for i = sigma1range
        for j = sigma2range
            idx = idx + 1;
        end
    end
    n = idx;

    reconstructions = zeros(msiN, n, 81);
    evaluation = zeros(msiN, n, 3);
    settings = struct('rho', [], 'windowDim', [], 'noiseType', {}, 'smoothingMatrixMethod', {}, 'noiseParam', {});
    w = warning('off', 'all');

    for k = 1:msiN

        % Retrieve MSI data

        infile = fullfile(getSetting('systemdir'), 'infiles', strcat('poi_', num2str(k), '.mat'));
        load(infile, 'poiName', 'poiRAW', 'poiSegmentMask', ...
            'roiSeeds', 'measuredSpectrum', 'poiWhite');

        msi = poiRAW;
        mask = poiSegmentMask;
        measured = measuredSpectrum;

        idx = 0;
        for rho = rhoRange
            for windowDim = windowRange
                for sigma = olympusSigmaRange
                    idx = idx + 1;
                    settings(idx).rho = rho;
                    settings(idx).windowDim = windowDim;
                    settings(idx).noiseType = 'spatiospectralolympus';
                    settings(idx).smoothingMatrixMethod = 'Cor_Sample';
                    settings(idx).noiseParam = {sigma};

                    setSetting('rho', rho);
                    setSetting('windowDim', windowDim);
                    setSetting('noiseType', 'spatiospectralolympus');
                    setSetting('smoothingMatrixMethod', 'Cor_Sample');
                    setSetting('noiseParam', sigma);

                    [estimated, gfc, nmse] = estimateReflectance(msi, mask, measured, ID(k));
                    reconstructions(k, idx, :) = estimated;
                    evaluation(k, idx, :) = [rmse, nmse, gfc];
                end
            end
        end


        for i = sigma1range
            idx = idx + 1;
            settings(idx).noiseType = 'sameForChannel';
            settings(idx).smoothingMatrixMethod = 'Cor_Sample';
            settings(idx).noiseParam = {i};

            setSetting('noiseType', 'sameForChannel');
            setSetting('smoothingMatrixMethod', 'Cor_Sample');
            setSetting('noiseParam', i);
            [estimated, gfc, nmse] = estimateReflectance(msi, mask, measured, ID(k));
            reconstructions(k, idx, :) = estimated;
            evaluation(k, idx, :) = [rmse, nmse, gfc];
        end

        for i = sigma1range
            for j = sigma2range
                idx = idx + 1;
                settings(idx).noiseType = 'spatial';
                settings(idx).smoothingMatrixMethod = 'Cor_Sample';
                settings(idx).noiseParam = {i, j};

                setSetting('noiseType', 'spatial');
                setSetting('smoothingMatrixMethod', 'Cor_Sample');
                setSetting('noiseParam', [i, j]);

                [estimated, gfc, nmse] = estimateReflectance(msi, mask, measured, ID(k));
                reconstructions(k, idx, :) = estimated;
                evaluation(k, idx, :) = [rmse, nmse, gfc];
            end
        end

    end

    warning(w);

    evaluationMeans = squeeze(mean(evaluation, 1));

    for i = 1:size(reconstructions, 2)
        settings(i).rmse = evaluationMeans(i, 1);
        settings(i).nmse = evaluationMeans(i, 2);
        settings(i).gfc = evaluationMeans(i, 3);
    end

    savefile = fullfile(getSetting('savedir'), 'Reconstruction Parameter Optimization', 'reconstructionComparison.mat');
    save(savefile, 'reconstructions', 'evaluation', 'settings', 'evaluationMeans');

else

    %% Load comparison results
    load('D:\temp\Google Drive\titech\research\output\saitama_v8_min_region_bright\Reconstruction Parameter Optimization\reconstructionComparison.mat');
end

setSetting('saveImages', true);

idx = find(cellfun(@(x) strcmp(x, 'sameForChannel'), {settings.noiseType}));
sigmas = cell2mat([settings(idx).noiseParam]);

setSettin('plotName', fullfile(getSetting('savedir'), 'Reconstruction Parameter Optimization', 'simpleWiener'));
plotNMSE(sigmas, evaluationMeans(idx, 2), '$$log_{10}( \sigma_1 )$$', ...
    'Average NMSE', 'NMSE for simple Wiener Reconstruction', 1);
setSettin('plotName', fullfile(getSetting('savedir'), 'Reconstruction Parameter Optimization', 'simpleWiener_gfc'));
plotGFC(sigmas, evaluationMeans(idx, 3), '$$log_{10}( \sigma_1 )$$', ...
    'Average GFC', 'GFC for simple Wiener Reconstruction', 2);

idx = find(cellfun(@(x) contains(x, 'spatiospectralolympus'), {settings.noiseType}));
rhos = [settings(idx).rho];
windowDim = [settings(idx).windowDim];
groups = findgroups(rhos, windowDim);
bestIdxs = zeros(1, max(groups));
for i = 1:max(groups)
    groupIdx = find(groups == i);
    [~, bestIdxs(i)] = min(evaluationMeans(idx(groupIdx), 2));
end
rhos = unique(rhos);
windowDim = unique(windowDim);

setSetting('plotName', fullfile(getSetting('savedir'), 'Reconstruction Parameter Optimization', 'Spatiospectral'));
plotNMSE2(rhos, windowDim, evaluationMeans(idx(bestIdxs), 2), 'Spatial Correlation', ...
    'Window Dimension', 'Average NMSE', 'NMSE for Spatiospectral Wiener Reconstruction', 1);
setSetting('plotName', fullfile(getSetting('savedir'), 'Reconstruction Parameter Optimization', 'Spatiospectral_gfc'));
plotGFC2(rhos, windowDim, evaluationMeans(bestIdxs, 3), 'Spatial Correlation', ...
    'Window Dimension', 'Average GFC', 'GFC for Spatiospectral Wiener Reconstruction', 2);

idx = find(cellfun(@(x) strcmp(x, 'spatial'), {settings.noiseType}));
sigmas = cellfun(@(x) (sqrt(0.5) * x{1} + x{2}), {settings(idx).noiseParam});
[sigmas, sorted] = sort(sigmas, 'ascend');
idx = idx(sorted);
setSetting('plotName', fullfile(getSetting('savedir'), 'Reconstruction Parameter Optimization', 'spatial'));
plotNMSE(sigmas, evaluationMeans(idx, 2), '$$log_{10}( \sigma_1 )$$', ...
    'Average NMSE', 'NMSE for simple Wiener Reconstruction', 3);
setSetting('plotName', fullfile(getSetting('savedir'), 'Reconstruction Parameter Optimization', 'spatial_gfc'));
plotGFC(sigmas, evaluationMeans(idx, 3), '$$log_{10}( \sigma_1 )$$', ...
    'Average GFC', 'GFC for simple Wiener Reconstruction', 4);

[~, idx] = sort([settings.nmse]);
best = settings(idx(1:10));

plotReconstructionGIF(Spectra, reconstructions);

bestSpatioSpect = 1;
bestSimple = 96;
setSetting('plotName', fullfile(getSetting('savedir'), 'Reconstruction Parameter Optimization', 'simpleHist'));
plotGFCHistogram(evaluation(:, bestSimple, 3), 1);
setSetting('plotName', fullfile(getSetting('savedir'), 'Reconstruction Parameter Optimization', 'spatiospectHist'));
plotGFCHistogram(evaluation(:, bestSpatioSpect, 3), 2);
