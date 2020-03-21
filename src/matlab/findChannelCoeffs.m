close all; clc; 
k = 20;
msiType = 'extended'; %'extended'; % 'max';
 
% for all pixels without bg removal 
removebg = false; 
[msi, ~, ~, ~, ~, ~] = getImage(k, options, msiType, removebg, true);
[coeff,latent,explained]  = doPixelPCA(msi);

% for only the pixels of the speciment 
removebg = true; 
[msi, ~, ~, ~, ~, ~] = getImage(k, options, msiType, removebg, true);
[coeff,latent,explained]  = doPixelPCA(msi);

normMsi = getNormMSI(msi);
[Loadings1,specVar1,T,stats] = factoran(normMsi,2,'rotate','none') %'Xtype','cov',


[msi, ~, ~, ~, ~, ~] = getImage(k, options, msiType, removebg, false);
plotMSI(msi, 1, options.saveOptions);
title('Msi (Extended)')
[difMsi, avgMsi] = getDifMSI(msi, 'toAverage') ;
plotMSI(avgMsi, 2, options.saveOptions);
title('Average')
plotMSI(difMsi, 3, options.saveOptions);
title('Difference to Average')
[difMsi, ~] = getDifMSI(msi, 'toNextBand') ;
plotMSI(difMsi, 4, options.saveOptions);
title('Difference to Next Band')
%xcorrMsi = getXcorrMSI(msi) ;
%plotMSI(xcorrMsi, 5, options.saveOptions);
%title('Difference to Next Band') 
[normMsi] = getNormMSI(msi);
plotMSI(normMsi, 6, options.saveOptions);
title('Normalized Msi') 

function [coeff,latent,explained]  = doPixelPCA(columns)

    [coeff,score,latent,tsquared,explained,mu] = pca(columns,'NumComponents',2, 'Centered', true);

    coeff
    latent
    explained

end 

function [difMsi, avgMsi] = getDifMSI(msi, type)
    if strcmp(type, 'toAverage')
        channels = size(msi, 1);
        avgMsi =  mean(msi,1);
        avgMsiRep = repmat(avgMsi, channels, 1, 1);   
        difMsi = bsxfun(@(x,y) x - y, msi, avgMsiRep); 
    elseif strcmp(type, 'toNextBand')
        [channels, m, n] = size(msi);
        difMsi = zeros(channels - 1, m, n);
        for i = 1:(channels - 1)
            difMsi(i, :, :) = squeeze(msi(i, :, :) - msi(i + 1, :, :));
        end 
        avgMsi = [];
    else 
        disp('Not supported method.');
    end    
end 

function [xcorrMsi] = getXcorrMSI(msi)
    [channels, m, n] = size(msi);
    for i = 1:(channels - 1)
        xcorrMsi(i, :, :) = xcorr2(squeeze(msi(i, :, :)), squeeze(msi(i + 1, :, :)));
    end
end

function [normMsi] = getNormMSI(msi)
    if ndims(msi) > 2 
        normMsi = zeros(size(msi));
        for i = 1:size(msi,1)
            normMsi(i, :, :) = squeeze(msi(i, :, :)) - mean(mean(squeeze(msi(i, :, :))));
        end 
    else 
        normMsi = zeros(size(msi));
        for i = 1:size(msi,2)
            normMsi(:,i) = squeeze(msi(:, i)) - mean(squeeze(msi(:, i)));
        end 
    end 
end

