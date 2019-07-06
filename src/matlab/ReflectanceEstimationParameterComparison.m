sigma1range = 10.^(-5:0.2:0);
sigma2range = 10.^(-5:0.2:0);
sigma1n = length(sigma1range);
sigma2n = length(sigma2range);
aRange = 0:0.2:1;
an = length(aRange);
gammaRange = 0.5:0.5:3;
gamman = length(gammaRange);

options.pixelValueSelectionMethod = 'extended';
n = 1051;
reconstructions = zeros(msiN, n, 81);
evaluation = zeros(msiN, n, 3);
settings = struct('rho', [], 'windowDim', [], 'noiseType', {}, 'smoothingMatrixMethod', {}, 'noiseParam', {});
%w = warning('off', 'all');

for k = 1:msiN
    
    % Retrieve MSI data
    
    infile = fullfile(options.systemdir, 'infiles', strcat('poi_', num2str(k), '.mat'));
    load(infile, 'poiName', 'poiRAW', 'poiSegmentMask', ...
        'roiSeeds', 'measuredSpectrum', 'poiWhite');
        
    msi = poiRAW;
    mask = poiSegmentMask;
    measured = measuredSpectrum;
    
    idx = 0;
    for rho = [0.6, 0.9, 0.95, 0.97, 0.985]
        for windowDim = [ 3, 5, 9 ]
            for sigma = kron( 10.^(-5:-3), [1, 5])
                idx = idx + 1;
                settings(idx).rho = rho;
                settings(idx).windowDim = windowDim;
                settings(idx).noiseType = 'spatiospectralolympus';
                settings(idx).smoothingMatrixMethod = 'Cor_Sample';
                settings(idx).noiseParam = {sigma};
                
                options.rho = rho;
                options.windowDim = windowDim; 
                options.noiseType = 'spatiospectralolympus';
                options.smoothingMatrixMethod = 'Cor_Sample';
                options.noiseParam = sigma;
     
                [estimated, rmse, nmse, ~] = reflectanceEstimation(msi, mask, measured, ID(k), options); 
                gfc = GoodnessOfFit(estimated, measured);
                reconstructions(k, idx, :) = estimated;
                evaluation(k, idx, :) = [rmse, nmse, gfc];
            end
        end
    end
    
    for i = 1:sigma1n
        idx = idx + 1;
        settings(idx).noiseType = 'sameForChannel';
        settings(idx).smoothingMatrixMethod = 'Cor_Sample';
        settings(idx).noiseParam = {sigma1range(i)};
                
        options.noiseType = 'sameForChannel';
        options.smoothingMatrixMethod = 'Cor_Sample';
        options.noiseParam = sigma1range(i);
        [estimated, rmse, nmse, ~] = reflectanceEstimation(msi, mask, measured, ID(k), options); 
        gfc = GoodnessOfFit(estimated, measured);
        reconstructions(k, idx, :) = estimated;
        evaluation(k, idx, :) = [rmse, nmse, gfc];
    end
    
    for i = 1:sigma1n
        for j = 1:sigma2n
            idx = idx + 1;
            settings(idx).noiseType = 'spatial';
            settings(idx).smoothingMatrixMethod = 'Cor_Sample';
            settings(idx).noiseParam = {sigma1range(i) , sigma2range(j)};
        
            options.noiseType = 'spatial';
            options.smoothingMatrixMethod = 'Cor_Sample';
            options.noiseParam = [sigma1range(i) , sigma2range(j)];
            [estimated, rmse, nmse, ~] = reflectanceEstimation(msi, mask, measured, ID(k), options); 
            gfc = GoodnessOfFit(estimated, measured);
            reconstructions(k, idx, :) = estimated;
            evaluation(k, idx, :) = [rmse, nmse, gfc];
        end
    end
    
    for i = -3:3
        idx = idx + 1;
        
        settings(idx).noiseType = 'fromOlympus';
        settings(idx).smoothingMatrixMethod = 'Cor_Sample';
        settings(idx).noiseParam = {10^i};
            
        options.noiseType = 'fromOlympus';
        options.noiseParam = 10^i;
        
        [estimated, rmse, nmse, ~] = reflectanceEstimation(msi, mask, measured, ID(k), options); 
        gfc = GoodnessOfFit(estimated, measured);
        reconstructions(k, idx, :) = estimated;
        evaluation(k, idx, :) = [rmse, nmse, gfc];
    end
    
    for a = aRange
        for gamma = gammaRange
            for i = -3:3              
                idx = idx + 1;
                settings(idx).noiseType = 'fromOlympus';
                settings(idx).smoothingMatrixMethod = 'adaptive';
                settings(idx).noiseParam = {a, gamma, 10^i};

                options.alpha = a;
                options.gamma = gamma;
                options.smoothingMatrixMethod = 'adaptive';
                options.noiseType = 'fromOlympus';
                options.noiseParam = 10^i;

                [estimated, rmse, nmse, ~] = reflectanceEstimation(msi, mask, measured, ID(k), options); 
                gfc = GoodnessOfFit(estimated, measured);
                reconstructions(k, idx, :) = estimated;
                evaluation(k, idx, :) = [rmse, nmse, gfc];
            end
        end
    end
    
end
warning(w);

evaluationMeans = squeeze(mean(evaluation, 1));

for i = 1:size(reconstructions, 2)
    settings(i).rmse = evaluationMeans(i,1);
    settings(i).nmse = evaluationMeans(i,2);
    settings(i).gfc = evaluationMeans(i,3);
end

savefile = fullfile(options.saveOptions.savedir, 'RecostructionComparison', 'reconstructionComparison.mat');
save(savefile, 'reconstructions', 'evaluation', 'settings', 'evaluationMeans');


saveOptions.saveImages = true;

idx = find( cellfun(@(x) strcmp(x, 'sameForChannel'), {settings.noiseType}));
sigmas = cell2mat([settings(idx).noiseParam]);

saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Reconstruction Parameter Optimization','simpleWiener');
plotNMSE(sigmas, evaluationMeans(idx, 2), '$$log_{10}( \sigma_1 )$$',...
    'Average NMSE', 'NMSE for simple Wiener Reconstruction', 1, saveOptions);
saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Reconstruction Parameter Optimization','simpleWiener_gfc');
plotGFC(sigmas, evaluationMeans(idx, 3), '$$log_{10}( \sigma_1 )$$',...
    'Average GFC', 'GFC for simple Wiener Reconstruction', 2, saveOptions);

idx = find( cellfun(@(x) strcmp(x, 'spatiospectralolympus'), {settings.noiseType}));
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

saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Reconstruction Parameter Optimization','Spatiospectral');
plotNMSE2(rhos, windowDim, evaluationMeans(idx(bestIdxs), 2), 'Spatial Correlation',...
    'Window Dimension', 'Average NMSE', 'NMSE for Spatiospectral Wiener Reconstruction', 1, saveOptions);
saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Reconstruction Parameter Optimization','Spatiospectral_gfc');
plotGFC2(rhos, windowDim, evaluationMeans(bestIdxs, 3),  'Spatial Correlation',...
    'Window Dimension','Average GFC', 'GFC for Spatiospectral Wiener Reconstruction', 2, saveOptions);

idx = find( cellfun(@(x) strcmp(x, 'spatial'), {settings.noiseType}));
sigmas = cellfun(@(x) (sqrt(0.5) * x{1} + x{2}), {settings(idx).noiseParam});
[sigmas, sorted] = sort(sigmas, 'ascend');
idx = idx(sorted);
saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Reconstruction Parameter Optimization','spatial');
plotNMSE(sigmas, evaluationMeans(idx, 2), '$$log_{10}( \sigma_1 )$$',...
    'Average NMSE', 'NMSE for simple Wiener Reconstruction', 1, saveOptions);
saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Reconstruction Parameter Optimization','spatial_gfc');
plotGFC(sigmas, evaluationMeans(idx, 3), '$$log_{10}( \sigma_1 )$$',...
    'Average GFC', 'GFC for simple Wiener Reconstruction', 2, saveOptions);

[~, idx] =  sort([settings.nmse]);
best = settings(idx(1:10));

plotReconstructionGIF(Spectra, reconstructions);
    
    
function gfc = GoodnessOfFit(reconstructed, measured)
    gfc = abs(reconstructed * measured') / ( sqrt(sum(reconstructed.^2)) * sqrt(sum(measured.^2)));
end
