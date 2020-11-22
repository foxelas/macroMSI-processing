
%% Preprocessing for Maps
setSetting('saveImages', false);
preprocessingForMaps;
% 
% %% Chromophore Absorbance Maps
resultFixed = prepareMaps(true, preprocFixed, fixedWhiteReference, fixedMask, ID(fixedId), fixedRoiCorners);
resultUnfixed = prepareMaps(false, preprocUnfixed, unfixedWhiteReference, unfixedMask, ID(unfixedId), unfixedRoiCorners);

% setSetting('saveImages', false);
%% Plot ROIs
figTitle = strjoin({tissueStates{1}, 'vs', tissueStates{2}, strcat('(', lesion, ')')}, ' ' );
setSetting('plotName', fullfile(getSetting('savedir'), getSetting('common'), strcat(lesion, '_', 'rois.png')));
plotFunWrapper(1, @plotMontage, resultUnfixed.Iaug, resultFixed.Iaug, figTitle);

%% Plot Maps Collectively
plotFunWrapper(2, @plotMontageScaled, makeImageList(resultUnfixed.MelMaps, resultFixed.MelMaps), ...
    {mapData.Mask}, imageNames, strcat(lesion, '_Full_Mel'));
plotFunWrapper(3, @plotMontageScaled, makeImageList(resultUnfixed.HbMaps, resultFixed.HbMaps), ...
    {mapData.Mask}, imageNames, strcat(lesion, '_Full_Hb'));

for k = 1:length(roiNames)
    masks = {ones(size(resultUnfixed.RoiMelMaps(1, k))), ones(size(resultUnfixed.RoiMelMaps(2, k)))};
    plotFunWrapper(4, @plotMontageScaled, makeImageList(resultUnfixed.RoiMelMaps(:, k), resultFixed.RoiMelMaps(:, k)), ...
        masks, imageNames, strcat(lesion, '_', roiNames{k}, '_Mel'));
    plotFunWrapper(5, @plotMontageScaled, makeImageList(resultUnfixed.RoiHbMaps(:, k), resultFixed.RoiHbMaps(:, k)), ...
        masks, imageNames, strcat(lesion, '_', roiNames{k}, '_Hb'));
end

%% Get Similarity Metrics
close all;
for i = 1:length(mapMethods)
    indexes = [mapData.Index];
    [metricsMel(i, :)] = getImageSimilarity(resultUnfixed.MelMaps{i}, resultFixed.MelMaps{i}, strcat(mapMethods{i}, '_mel'), lesion, tissueStates);
    [metricsHb(i, :)] = getImageSimilarity(resultUnfixed.HbMaps{i}, resultFixed.HbMaps{i}, strcat(mapMethods{i}, '_hb'), lesion, tissueStates);

    for k = 1:length(roiNames)

        roiMetricsMel(i, k, :) = getImageSimilarity(resultUnfixed.RoiMelMaps{i, k}, resultFixed.RoiMelMaps{i, k}, ...
            strcat(mapMethods{i}, '_', roiNames{k}, '_mel'), lesion, tissueStates);
        roiMetricsHb(i, k, :) = getImageSimilarity(resultUnfixed.RoiHbMaps{i, k}, resultFixed.RoiHbMaps{i, k}, ...
            strcat(mapMethods{i}, '_', roiNames{k}, '_hb'), lesion, tissueStates);
    end

end

metricsMelTable = makeTable(metricsMel, metricsNames, mapMethods');
metricsHbTable = makeTable(metricsHb, metricsNames, mapMethods');

roiMetricsMelTable = makeTable(roiMetricsMel, metricsNames, mapMethods', roiNames');
roiMetricsHbTable = makeTable(roiMetricsHb, metricsNames, mapMethods', roiNames');

%% Export latex tables 
metricsMelTableText = table2latex(metricsMelTable, [1, 3:6], 'mel', strjoin({'Similarity of Melanin Maps for', lesion}, ' '));
metricsHbTableText = table2latex(metricsHbTable, [1, 3:6], 'hb', strjoin({'Similarity of Hemoglobin Maps for', lesion}, ' '));
roiMetricsMelTableText = table2latex(roiMetricsMelTable, [1:2, 4:7], 'melroi', strjoin({'Similarity of Melanin Maps for', lesion , 'ROIs'}, ' '));
roiMetricsHbTableText = table2latex(roiMetricsHbTable, [1:2, 4:7], 'hbroi', strjoin({'Similarity of Hemoglobin Maps for', lesion, 'ROIs'}, ' '));


%% Functions
function T = makeTable(metrics, names, methods, rois)
if nargin < 4
    T = array2table(metrics);
    T.Properties.VariableNames(1:length(names)) = names;
    T.Method = methods;
    T = [T(:, end), T(:, 1:end-1)];
else
    metrics = reshape(metrics, [size(metrics, 1) * size(metrics, 2), size(metrics, 3)]);
    T = array2table(metrics);
    T.Properties.VariableNames(1:length(names)) = names;
    T.Method = repmat(methods, [size(metrics, 1) / numel(methods), 1]);
    T.RoiType = reshape(repmat(rois, [1, numel(methods)])', [numel(methods) * numel(rois), 1]);
    T = [T(:, (end -1):end), T(:, 1:end-2)];
end

end

function imageList = makeImageList(images1, images2)
imageList = cell(numel(images1)*2, 1);
for i = 1:(numel(images1))
    imageList{i*2-1} = images1{i};
    imageList{i*2} = images2{i};
end
end

