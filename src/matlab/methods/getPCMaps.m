function [pcComps] = getPCMaps(id, msiType, normType, removeBg, tform2, newDims)

nComps = 3;
savedir = getSetting('savedir');
mapdir = getSetting('pca');

% [msi, ~, specimenMask, height, width, channels] = getImage(k, msiType, removeBg, false, normType);
[colMsi, ~, colSpecimenMask, height, width, channels] = getImage(id, msiType, removeBg, true, normType, tform2, newDims);
specimenIdxs = find(colSpecimenMask == 1)';
specimenMask = reshape(colSpecimenMask, height, width);


%[coeff]  = doPixelPCA(colMsi);
[coeff] = applyOnQualityPixels(@doPixelPCA, colMsi);

centeredX = colMsi - mean(colMsi);
sc = centeredX * coeff;
pcComps = zeros(nComps, height, width);
for i = 1:nComps
    pcComp = zeros(length(colSpecimenMask), 1);
    pcComp(specimenIdxs) = sc(:, i);
    pcCompImage = reshape(pcComp, height, width);
    %     figure(i);
    %     imagesc(pcCompImage);
    %     colorbar;
    %     title(sprintf('Principal Component %d', i))
    %     axis off;
    %     %caxis([-2,2])

    pcComps(i, :, :) = pcCompImage;
    setSetting('plotName', fullfile(savedir, mapdir, strcat('norm_', normType), strcat(num2str(id), '_PC_', num2str(i), '.png')));
    plots(i, @plotMap, pcCompImage, specimenMask, [], false, strcat('PC', num2str(i)));
end

end


function [coeff, scores, latent, explained] = doPixelPCA(columns)

[coeff, scores, latent, tsquared, explained, mu] = pca(columns, 'NumComponents', 4, 'Centered', true);

coeff
latent'
explained'

end