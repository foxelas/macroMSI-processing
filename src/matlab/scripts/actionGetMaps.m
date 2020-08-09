% normType = 'divMacbeth';
% removeBg = 'true';
% fixedId = 27; 
% unfixedId = 28; 
%  
% % %% Preprocessing for PCA 
% % 
% % % msiType = 'extended'; 
% % % preprocessingForMaps; 
% % % 
% % % %% PCA maps 
% % % 
% % % close all; pcCompsFixed = getPCMaps(fixedId, msiType, normType, removeBg, tform2, newDims);
% % % close all; pcCompsUnfixed = getPCMaps(unfixedId, msiType, normType, removeBg, tform2, newDims);
% % % 
% % % %% PCA Similarity 
% % % close all;
% % % for i = 1: 3
% % %     pc_fixed = squeeze(pcCompsFixed(i,:,:)); 
% % %     pc_unfixed = squeeze(pcCompsUnfixed(i,:,:)); 
% % %     [ssimval, cc] = getImageSimilarity(pc_unfixed, pc_fixed, strcat('_pc', ...
% % %         num2str(i)), getSetting('pca'), unfixedId, fixedId, name1, name2);
% % %     pause(0.5);
% % % end 
% 
% %% Preprocessing for Maps  
% 
% msiType = 'adjusted'; 
% preprocessingForMaps; 

%% Chromophore Absorbance Maps 

metricsNames = {'SSIM', 'SSD', 'NCC', 'CC', 'CR', 'SimilarityHistIntersection', 'KLDivergece', 'EarthMoversDistance'}; 
mapMethods = { 'vasefi', 'ding', 'ours', 'diebele', 'kapsokalyvas', 'kuzmina'};
metricsMel = zeros(length(mapMethods), length(metricsNames));
metricsHb = zeros(length(mapMethods), length(metricsNames));

roiNames = {'Hb', 'Norm', 'Mel'};
roiCorners = {[316, 382, 242, 295], [159, 252, 167, 214], [398, 440, 83, 137]};
roiColors = {'m', 'g', 'c'}; 
% cornersHb = [316, 382, 242, 295];   %bboxHb = [316 242 382-316 295-242];
% cornersNorm = [159, 252, 167, 214];    %bboxNorm = [159 154 252-159 225-154];
% cornersMel = [398, 440, 83, 137];   %bboxMel = [398 83 440-398 137-83];

masks = {unfixedMask, fixedMask};
msis = {preprocUnfixed, preprocFixed};
whites = {unfixedWhiteReference, fixedWhiteReference};
mapData = struct('Msi', [], 'Mask', [], 'White', [], 'Iaug', [], 'RoiMasks', [], 'RoiBboxes', []);
for j = 1:2
    mapData(j).Msi = msis{j};
    mapData(j).Mask = masks{j};
    mapData(j).White = whites{j};
    Iaug = mapData(j).White;
    roiMask = cell(numel(roiNames), 1);
    bbox = cell(numel(roiNames), 1);
    for i = 1:numel(roiNames)
        [roiMask{i}, bbox{i}] = getBoundingBoxMask(roiCorners{i},  masks{j});
        Iaug = insertShape(Iaug,'rectangle', bboxMel,'LineWidth',5, 'Color', roiColors{i});
    end
    mapData(j).RoiMasks = roiMask;
    mapData(j).RoiBboxes = bbox;
    mapData(j).Iaug = Iaug;
    figure(j); imshow(Iaug); savePlot(j);
end
% [maskHb, bboxHb] = getBoundingBoxMask(cornersHb, mask);
% [maskNorm, bboxNorm] = getBoundingBoxMask(cornersNorm, mask);
% [maskMel, bboxMel] = getBoundingBoxMask(cornersMel, mask);
% Iaug = insertShape(I,'rectangle', bboxMel,'LineWidth',5, 'Color', 'c');
% Iaug = insertShape(Iaug,'rectangle', bboxHb,'LineWidth',5, 'Color', 'm');
% Iaug = insertShape(Iaug,'rectangle', bboxNorm,'LineWidth',5, 'Color', 'g');
% imshow(Iaug)

setSetting('saveImages', false);

for i = 1:length(mapMethods)
    
    msi = preprocUnfixed;
    mask = unfixedMask;
    close all; [melMap1, hbMap1] = getMap(msi, mapMethods{i}, mask, ID(unfixedId));
    
    msi = preprocFixed;
    mask = fixedMask;
    close all; [melMap2, hbMap2] = getMap(msi, mapMethods{i}, mask, ID(fixedId));
    
    [metricsMel(i,:)] = getImageSimilarity(melMap1, melMap2, strcat(mapMethods{i}, '_mel'), unfixedId, fixedId, 'unfixed', 'fixed');
    [metricsHb(i,:)] = getImageSimilarity(hbMap1, hbMap2, strcat(mapMethods{i}, '_hb'), unfixedId, fixedId, 'unfixed', 'fixed');
    

    metrics = getImageSimilarity(imcrop(hbMap1, bboxHb), imcrop(hbMap2, bboxHb), [], unfixedId, fixedId, 'unfixed', 'fixed');

end 
metricsMelTable = makeTable(metricsMel, metricsNames, mapMethods');
metricsHbTable = makeTable(metricsHb, metricsNames, mapMethods');

function T = makeTable(metrics, names, methods)
    T = array2table(metrics);
    T.Properties.VariableNames(1:length(names)) = names;
    T.MappingMethod = methods;
    T = [T(:,end) T(:,1:end-1)];
end 


