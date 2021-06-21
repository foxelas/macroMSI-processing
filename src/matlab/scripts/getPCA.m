
%% Preprocessing for PCA

msiType = 'extended';
preprocessingForMaps;

%% PCA maps

close all; pcCompsFixed = getPCMaps(fixedId, msiType, normType, removeBg, tform2, newDims);
close all; pcCompsUnfixed = getPCMaps(unfixedId, msiType, normType, removeBg, tform2, newDims);

%% PCA Similarity
close all;
for i = 1:3
    pc_fixed = squeeze(pcCompsFixed(i, :, :));
    pc_unfixed = squeeze(pcCompsUnfixed(i, :, :));
    [ssimval, cc] = getImageSimilarity(pc_unfixed, pc_fixed, strcat('pc', ...
        num2str(i)), lesion, tissueStates);
    pause(0.5);
end
