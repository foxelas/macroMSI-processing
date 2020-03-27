close all; clc; 
k = 20;
msiType = 'extended'; %'extended'; % 'max';
 
% % for all pixels without bg removal 
% removebg = false; 
% [msi, ~, ~, ~, ~, ~] = getImage(k, options, msiType, removebg, true);
% [coeff, scores, latent,explained]  = doPixelPCA(msi);

% % for only the pixels of the specimen
% [msi, ~, specimenMask, height, width, ~] = getImage(k, options, msiType, false, true);
specimenIdxs = find(specimenMask == 1)';
colMsi = msi(specimenIdxs, :);
[coeff,scores, latent,explained]  = doPixelPCA(colMsi);
% [coeff,scores, latent,explained] = applyOnQualityPixels( @doPixelPCA, colMsi);

for i = 1 : 3
    pcComp = zeros(size(msi,1), 1);
    pcComp(specimenIdxs) = scores(:,i);
    pcCompImage = reshape(pcComp, height, width);
    figure(i);
    imagesc(pcCompImage);
    colorbar;
    title(sprintf('Principal Component %d', i))
    %caxis([-2,2])
end 

% pc1 healthy skin, pc2 hemoglobin, pc3 


% normMsi = getNormMSI(msi);
% [Loadings1,specVar1,T,stats] = factoran(normMsi,2,'rotate','none') %'Xtype','cov',


function [coeff,scores, latent,explained]  = doPixelPCA(columns)

    [coeff,scores,latent,tsquared,explained,mu] = pca(columns,'NumComponents',5, 'Centered', true);

    coeff
    latent
    explained

end 


