sigma1range = 10.^(-5:0.2:0);
sigma2range = 10.^(-5:0.2:0);
sigma1n = length(sigma1range);
sigma2n = length(sigma2range);
aRange = 0:0.2:1;
an = length(aRange);
gammaRange = 0.5:0.5:3;
gamman = length(gammaRange);
clear rmses nmses gfcs rr ww ss
% rmses = zeros(msiN, 7*3*5); %zeros(msiN, sigma1n + sigma1n * sigma2n + 7 + an * gamman);
% nmses = zeros(msiN, 7*3*5); %zeros(msiN, sigma1n + sigma1n * sigma2n  + 7 + an * gamman);
% gfcs =  zeros(msiN, 7*3*5); %zeros(msiN, sigma1n + sigma1n * sigma2n  + 7 + an * gamman);

options.pixelValueSelectionMethod = 'extended';
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
                options.rho = rho;
                options.windowDim = windowDim; 
                options.noiseType = 'spatiospectralolympus';
                options.smoothingMatrixMethod = 'Cor_Sample';
                options.noiseParam = sigma;
                idx = idx + 1;
                [estimated, rmses(k,idx), nmses(k,idx), ~] = reflectanceEstimation(msi, mask, measured, ID(k), options); 
                gfcs(k,idx) = GoodnessOfFit(estimated, measured);
                rr(idx) = rho;
                ss(idx) = sigma;
                ww(idx) = windowDim;
            end
        end
    end
%     
%     for i = 1:sigma1n
%         options.noiseType = 'sameForChannel';
%         options.smoothingMatrixMethod = 'Cor_Sample';
%         options.noiseParam = sigma1range(i);
%         [estimated, rmses(k,i), nmses(k,i), ~] = reflectanceEstimation(msi, mask, measured, ID(k), options);   
%         gfcs(k,i) = GoodnessOfFit(estimated, measured);
%     end
%     
%     idx = i;
%     for i = 1:sigma1n
%         for j = 1:sigma2n
%             idx = idx + 1;
%             options.noiseType = 'spatial';
%             options.smoothingMatrixMethod = 'Cor_Sample';
%             options.noiseParam = [sigma1range(i) , sigma2range(j)];
%             [estimated, rmses(k,idx), nmses(k,idx), ~] = reflectanceEstimation(msi, mask, measured, ID(k), options);
%             gfcs(k,idx) = GoodnessOfFit(estimated, measured);
%         end
%     end
%     
%     for i = -3:3
%         idx = idx + 1;
%         options.noiseType = 'fromOlympus';
%         options.noiseParam = 10^i;
%         [estimated, rmses(k,idx), nmses(k,idx), ~] = reflectanceEstimation(msi, mask, measured, ID(k), options);  
%         gfcs(k,idx) = GoodnessOfFit(estimated, measured);
%     end
%     
%     for a = aRange
%         for gamma = gammaRange
%             idx = idx + 1;
%             options.alpha = a;
%             options.gamma = gamma;
%             options.smoothingMatrixMethod = 'adaptive';
%             options.noiseType = 'fromOlympus';
%             options.noiseParam = 0.1;
%             [estimated, rmses(k,idx), nmses(k,idx), ~] = reflectanceEstimation(msi, mask, measured, ID(k), options);  
%             gfcs(k,idx) = GoodnessOfFit(estimated, measured);
%         end
%     end
%     
end
warning(w);

%save('rmsenmse.mat', 'rmses', 'nmses', 'gfcs');

saveOptions.saveImages = true;

% saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Reconstruction Parameter Optimization','simpleWiener');
% plotNMSE(sigma1range, mean(nmses(:, 1:sigma1n)), '$$log_{10}( \sigma_1 )$$',...
%     'Average NMSE', 'Error for simple Wiener Reconstruction', 1, saveOptions);
% saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Reconstruction Parameter Optimization','simpleWiener_gfc');
% plotGFC(sigma1range, mean(gfcs(:, 1:sigma1n)), '$$log_{10}( \sigma_1 )$$',...
%     'Average GFC', 'GFC for simple Wiener Reconstruction', 2, saveOptions);
% 
% saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Reconstruction Parameter Optimization','spatialDenoising');
% kk = (sigma1n + 1):(sigma1n * sigma2n + sigma1n);
% plotNMSE2(sigma1range, sigma2range, mean(nmses(:,kk)), '$$log_{10}( \sigma_1 )$$',...
%     '$$log_{10}( \sigma_2 )$$', 'Average NMSE', 'Error for Spatially Denoised Wiener Reconstruction', 3, saveOptions);
% 
% saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Reconstruction Parameter Optimization','spatialDenoising_gfc');
% kk = (sigma1n + 1):(sigma1n * sigma2n + sigma1n);
% plotGFC2(sigma1range, sigma2range, mean(gfcs(:,kk)), '$$log_{10}( \sigma_1 )$$',...
%     '$$log_{10}( \sigma_2 )$$', 'Average GFC', 'GFC for Spatially Denoised Wiener Reconstruction', 4, saveOptions);
% 
% ii = 0;
% xx = zeros(1, sigma1n);
% for i =1:sigma1n
%     for j=1:sigma2n
%         ii = ii + 1;
%         xx(ii) = sqrt(0.5) * sigma1range(i) + sigma2range(j);
%     end
% end
% [xx, idxx] = sort(xx, 'ascend');
% y = mean(gfcs(:,kk));
% yy = y(:,idxx);
% saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Reconstruction Parameter Optimization','spatialDenoising_gfc_flat');
% plotGFC(xx, yy, '$$log_{10}(\sqrt{0.5} \sigma_1  + \sigma_2)$$',...
%     'Average GFC', 'GFC for spatially Denoised Wiener Reconstruction', 5, saveOptions);
% y = mean(nmses(:,kk));
% yy = y(:,idxx);
% saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Reconstruction Parameter Optimization','spatialDenoising_flat');
% plotNMSE(xx, yy, '$$log_{10}(\sqrt{0.5} \sigma_1  + \sigma_2)$$',...
%     'Average NMSE', 'Error for spatially Denoised Wiener Reconstruction', 6, saveOptions);
% 
% saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Reconstruction Parameter Optimization','olympus_gfc');
% plotGFC(-3:3,  mean(gfcs(:,(sigma1n * sigma2n + sigma1n + 1):(sigma1n * sigma2n + sigma1n + 7))), '$$10^{n}$$ noise multiplier',...
%     'Average GFC', 'GFC for Wiener Reconstruction with noise model from Olympus', 7, saveOptions);
% 
% saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Reconstruction Parameter Optimization','adaptive_gfc');
% adaptStartIdx = sigma1n * sigma2n + sigma1n + 7 + 1;
% plotGFC2(aRange, gammaRange, mean(gfcs(:,adaptStartIdx:end)), 'alpha',...
%     'gamma', 'Average GFC', 'GFC for Adaptive Wiener Reconstruction', 8, saveOptions);
% saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Reconstruction Parameter Optimization','adaptive');
% adaptStartIdx = sigma1n * sigma2n + sigma1n + 7 + 1;
% plotNMSE2(aRange, gammaRange, mean(nmses(:,adaptStartIdx:end)), 'alpha',...
%     'gamma', 'Average NMSE', 'Error for Adaptive Wiener Reconstruction', 9, saveOptions);
% 

%     idx = 0;
%     for rho = [0.6, 0.9, 0.95, 0.97, 0.985]
%         for windowDim = [3, 5, 9]
%             for sigma = 10.^(-3:0.5:0)
%                 idx = idx +1;
%                 rr(idx) = rho;
%                 ww(idx) = windowDim;
%                 sw(idx) = sigma;
%             end
%         end
%     end
    
    
function gfc = GoodnessOfFit(reconstructed, measured)
    gfc = abs(reconstructed * measured') / ( sqrt(sum(reconstructed.^2)) * sqrt(sum(measured.^2)));
end
