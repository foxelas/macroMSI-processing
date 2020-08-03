normType = 'divMacbeth';
removeBg = 'true';
fixedId = 27; 
unfixedId = 28; 
 
% %% Preprocessing for PCA 
% 
% % msiType = 'extended'; 
% % preprocessingForMaps; 
% % 
% % %% PCA maps 
% % 
% % close all; pcCompsFixed = getPCMaps(fixedId, msiType, normType, removeBg, tform2, newDims);
% % close all; pcCompsUnfixed = getPCMaps(unfixedId, msiType, normType, removeBg, tform2, newDims);
% % 
% % %% PCA Similarity 
% % close all;
% % for i = 1: 3
% %     pc_fixed = squeeze(pcCompsFixed(i,:,:)); 
% %     pc_unfixed = squeeze(pcCompsUnfixed(i,:,:)); 
% %     [ssimval, cc] = getImageSimilarity(pc_unfixed, pc_fixed, strcat('_pc', ...
% %         num2str(i)), getSetting('pca'), unfixedId, fixedId, name1, name2);
% %     pause(0.5);
% % end 

%% Preprocessing for Maps  

msiType = 'adjusted'; 
preprocessingForMaps; 

%% Chromophore Absorbance Maps 

metricsNames = {'SSIM', 'SSD', 'NCC', 'CC', 'CR', 'SimilarityHistIntersection', 'KLDivergece', 'EarthMoversDistance'}; 
mapMethods = { 'vasefi', 'ding', 'ours', 'diebele', 'kapsokalyvas', 'kuzmina'};
metricsMel = zeros(length(mapMethods), length(metricsNames));
metricsHb = zeros(length(mapMethods), length(metricsNames));

for i = 1:length(mapMethods)
    
    msi = preprocUnfixed;
    mask = unfixedMask;
    close all; [melMap1, hbMap1] = getMap(msi, mapMethods{i}, mask, ID(unfixedId));
    
    msi = preprocFixed;
    mask = fixedMask;
    close all; [melMap2, hbMap2] = getMap(msi, mapMethods{i}, mask, ID(fixedId));
    
    [metricsMel(i,:), ~] = getImageSimilarity(melMap1, melMap2, strcat(mapMethods{i}, '_mel'), unfixedId, fixedId, 'unfixed', 'fixed');
    [metricsHb(i,:), ~] = getImageSimilarity(hbMap1, hbMap2, strcat(mapMethods{i}, '_hb'), unfixedId, fixedId, 'unfixed', 'fixed');
end 
metricsMelTable = makeTable(metricsMel, metricsNames, mapMethods');
metricsHbTable = makeTable(metricsHb, metricsNames, mapMethods');

function T = makeTable(metrics, names, methods)
    T = array2table(metrics);
    T.Properties.VariableNames(1:length(names)) = names;
    T.MappingMethod = methods;
    T = [T(:,end) T(:,1:end-1)];
end 


