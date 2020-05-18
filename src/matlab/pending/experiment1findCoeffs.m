close all; clc; 
k =  17; %20;
msiType = 'extended'; %'extended'; % 'max';
options.saveOptions.saveImages = false; 

% % for all pixels without bg removal 
% removebg = false; 
% [msi, ~, ~, ~, ~, ~] = getImage(k, options, msiType, removebg, true);
% [coeff, scores, latent,explained]  = doPixelPCA(msi);

[msiFull] = getImage(k, options, msiType, true);
% plotMSI(msiFull, 1, options.saveOptions);
% 
% for only the pixels of the specimen
[msi, ~, specimenMask, height, width, channels] = getImage(k, options, msiType, false, true);
[msi, ~] = getDifMSI(msi, 'toAverage') ;

specimenIdxs = find(specimenMask == 1)';
colMsi = msi(specimenIdxs, :);
%[coeff]  = doPixelPCA(colMsi);
[coeff] = applyOnQualityPixels( @doPixelPCA, colMsi);

centeredX = colMsi - mean(colMsi);
sc = centeredX * coeff; 

for i = 1 : 3
    pcComp = zeros(size(msi,1), 1);
    pcComp(specimenIdxs) = sc(:,i);
    pcCompImage = reshape(pcComp, height, width);
    figure(i);
    imagesc(pcCompImage);
    colorbar;
    title(sprintf('Principal Component %d', i))
    axis off;

    %caxis([-2,2])
end 

kCoeff = 2;
for i = 1 : 9
    w = zeros(channels,channels);
    w(i,i) = kCoeff; 
    enhancedImage = zeros(size(msi,1), channels);
    healthySkinComponent = sum(sc(:,1:2),2) + mean(colMsi);
    enhancedImage(specimenIdxs, :)  = (colMsi - healthySkinComponent) * w + colMsi ;
    enhancedImage = reshape(enhancedImage', channels, height, width);
    plotMSI(enhancedImage, 2, options.saveOptions);
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


