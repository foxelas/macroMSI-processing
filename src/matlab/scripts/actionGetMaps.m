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

roiNames = {'Hb', 'Norm', 'Mel'};
roiCorners = {[316, 382, 242, 295], [159, 252, 167, 214], [398, 440, 83, 137]};
roiColors = {'m', 'g', 'c'}; 
mapMethods =  { 'vasefi', 'ding', 'diebele', 'kapsokalyvas'};
 %{ 'vasefi', 'ding', 'ours', 'diebele', 'kapsokalyvas', 'kuzmina'};

masks = {unfixedMask, fixedMask};
msis = {preprocUnfixed, preprocFixed};
whites = {unfixedWhiteReference, fixedWhiteReference};
indexes = [unfixedId, fixedId];
imageNames = {'unfixed', 'fixed'};
mapData = struct('Name', [], 'Msi', [], 'Mask', [], 'White', [], 'Iaug', [], 'RoiMasks', [], 'RoiBboxes', [], 'ID', [], 'Index', [], 'MelMaps', [], 'HbMaps', [], 'RoiMelMaps', [], 'RoiHbMaps', []);
for j = 1:2
    mapData(j).Name = imageNames{j};
    mapData(j).Msi = msis{j};
    mapData(j).Mask = masks{j};
    mapData(j).White = whites{j};
    mapData(j).Index = indexes(j);
    mapData(j).ID = ID(indexes(j));
    
    Iaug = mapData(j).White;
    roiMask = cell(numel(roiNames), 1);
    bbox = cell(numel(roiNames), 1);
    for i = 1:numel(roiNames)
        [roiMask{i}, bbox{i}] = getBoundingBoxMask(roiCorners{i},  masks{j});
        Iaug = insertShape(Iaug,'rectangle', bbox{i},'LineWidth',5, 'Color', roiColors{i});
    end
    mapData(j).RoiMasks = roiMask;
    mapData(j).RoiBboxes = bbox;
    mapData(j).Iaug = Iaug;
    
    melMaps = cell(numel(mapMethods), 1);
    hbMaps = cell(numel(mapMethods), 1);
    roiMelMaps = cell(numel(mapMethods), numel(roiNames));
    roiHbMaps = cell(numel(mapMethods), numel(roiNames));
    for i = 1:length(mapMethods)
        msi =  msis{j};
        mask = masks{j};
        id =  ID(indexes(j));
        close all; [melMaps{i}, hbMaps{i}] = getMap(msi, mapMethods{i}, mask, id);     
        for k = 1:numel(roiNames)
            roiMelMaps{i, k} = imcrop(melMaps{i}, mapData(j).RoiBboxes{k});
            roiHbMaps{i, k} = imcrop(hbMaps{i}, mapData(j).RoiBboxes{k});
        end
        %% GetMap quality 
        ;
            
    end 
    mapData(j).MelMaps = melMaps;
    mapData(j).HbMaps = hbMaps;
    mapData(j).RoiMelMaps = roiMelMaps;
    mapData(j).RoiHbMaps = roiHbMaps;

end

% setSetting('saveImages', false);
metricsNames = {'SSIM', 'SSD', 'NCC', 'CC', 'CR', 'SimilarityHistIntersection', 'KLDivergece', 'EarthMoversDistance'}; 
metricsMel = zeros(length(mapMethods), length(metricsNames));
metricsHb = zeros(length(mapMethods), length(metricsNames));
roiMetricsMel = zeros(length(mapMethods), length(roiNames), length(metricsNames));
roiMetricsHb = zeros(length(mapMethods), length(roiNames), length(metricsNames));
imageNames =  makeNameList(mapMethods, {mapData.Name}); 
vsName = strcat(num2str(indexes(1)), 'vs', num2str(indexes(2)));

setSetting('saveImages', true);

%% Plot ROIs

setSetting('plotName', fullfile(getSetting('savedir'), getSetting('common'), strcat(vsName, '_', 'rois.png')));
plotFunWrapper(1, @plotMontage, mapData(1).Iaug, mapData(2).Iaug, strcat(mapData(1).Name, ' vs ', mapData(2).Name));

%% Plot Maps Collectively

plotFunWrapper(2, @plotMontageScaled, makeImageList(mapData(1).MelMaps, mapData(2).MelMaps), ...
    {mapData.Mask}, imageNames, strcat(vsName, '_Full_Mel'));
plotFunWrapper(3, @plotMontageScaled, makeImageList(mapData(1).HbMaps, mapData(2).HbMaps), ...
    {mapData.Mask}, imageNames, strcat(vsName, '_Full_Hb'));
    
for k = 1:length(roiNames)
    masks = {ones(size(mapData(1).RoiMelMaps(1,k))) ones(size(mapData(1).RoiMelMaps(2,k)))};
    plotFunWrapper(4, @plotMontageScaled, makeImageList(mapData(1).RoiMelMaps(:,k), mapData(2).RoiMelMaps(:,k)), ...
        masks, imageNames, strcat(vsName,'_', roiNames{k}, '_Mel'));
    plotFunWrapper(5, @plotMontageScaled, makeImageList(mapData(1).RoiHbMaps(:,k), mapData(2).RoiHbMaps(:,k)), ...
        masks, imageNames, strcat(vsName, '_', roiNames{k}, '_Hb'));  
end
  
%% Get Similarity Metrics 
close all; 
for i = 1:length(mapMethods)
    indexes = [mapData.Index];
    [metricsMel(i,:)] = getImageSimilarity(mapData(1).MelMaps{i}, mapData(2).MelMaps{i}, strcat(mapMethods{i}, '_mel'), indexes(1), indexes(2), mapData(1).Name, mapData(2).Name);
    [metricsHb(i,:)] = getImageSimilarity(mapData(1).HbMaps{i}, mapData(2).HbMaps{i}, strcat(mapMethods{i}, '_hb'), indexes(1), indexes(2), mapData(1).Name, mapData(2).Name);
    
    for k = 1:length(roiNames)
        
        roiMetricsMel(i,k,:) = getImageSimilarity(mapData(1).RoiMelMaps{i,k}, mapData(2).RoiMelMaps{i,k},...
            strcat(mapMethods{i}, '_', roiNames{k} , '_mel'), indexes(1), indexes(2), mapData(1).Name, mapData(2).Name);
        roiMetricsHb(i,k,:) = getImageSimilarity(mapData(1).RoiHbMaps{i,k}, mapData(2).RoiHbMaps{i,k},...
            strcat(mapMethods{i}, '_', roiNames{k} , '_hb'), indexes(1), indexes(2), mapData(1).Name, mapData(2).Name); 
    end
    
end 

metricsMelTable = makeTable(metricsMel, metricsNames, mapMethods');
metricsHbTable = makeTable(metricsHb, metricsNames, mapMethods');

roiMetricsHbTable = makeTable(roiMetricsHb, metricsNames, mapMethods', roiNames');
roiMetricsMelTable = makeTable(roiMetricsMel, metricsNames, mapMethods', roiNames');


%% Functions 
function T = makeTable(metrics, names, methods, rois)
    if nargin < 4 
        T = array2table(metrics);
        T.Properties.VariableNames(1:length(names)) = names;
        T.MappingMethod = methods;
        T = [T(:,end) T(:,1:end-1)];
    else
        metrics = reshape(metrics, [size(metrics,1) * size(metrics, 2), size(metrics, 3)]);
        T = array2table(metrics);
        T.Properties.VariableNames(1:length(names)) = names;
        T.MappingMethod = repmat(methods, [size(metrics,1) / numel(methods), 1]);
        T.RoiType = reshape(repmat(rois, [1 , numel(methods)])', [numel(methods)*numel(rois), 1]);
        T = [T(:,(end-1):end) T(:,1:end-2)];
    end
end 

function imageList = makeImageList(images1, images2)
    imageList = cell(numel(images1) * 2, 1);
    for i = 1:(numel(images1))
        imageList{i*2-1} = images1{i};
        imageList{i*2} = images2{i};
    end
end

function nameList = makeNameList(namesRows, namesColumns)
nameList = cell(numel(namesRows) * numel(namesColumns), 1);
for i = 1:numel(namesRows)
    for j = 1:numel(namesColumns)
        nameList{(i-1)*numel(namesColumns) + j} = strcat(namesColumns{j}, '_', namesRows{i});
    end
end

end 

