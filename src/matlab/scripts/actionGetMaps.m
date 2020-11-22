%% Initializing 
main; 

%% Preprocessing for Maps
fixedId = 13;
unfixedId = 14;
setSetting('saveImages', true);

%% Required settings
roiNames = getSetting('roiNames');
mapMethods = getSetting('mapMethods'); %{ 'vasefi', 'ding', 'ours', 'diebele', 'kapsokalyvas', 'kuzmina'};
metricsNames = getSetting('metricsNames');
tissueStates = getSetting('tissueStates');
normType = getSetting('normType');
removeBg = getSetting('removeBg');
msiType = getSetting('msiType');

roiCornerVals = delimread(fullfile(getSetting('systemdir'), getSetting('roiCornerFileName')), ',', {'num', 'text'});
textVals = roiCornerVals.text; 
roiCornerVals = roiCornerVals.num;
[unfixedRoiCorners, ~, ~] = getRoiCorners(unfixedId, roiCornerVals, textVals); 
[fixedRoiCorners, lesion, registrationType] = getRoiCorners(fixedId, roiCornerVals, textVals); 

%% Handle Fixed Image
[fixedMsi, fixedWhiteReference, fixedMask] = readAndNormalize(fixedId, ID, msiType, removeBg, normType);

%% Handle Unfixed Image
[unfixedMsi, unfixedWhiteReference, unfixedMask] = readAndNormalize(unfixedId, ID, msiType, removeBg, normType);

%% Register

if (~exist('tform', 'var'))
    if strcmp(registrationType, 'surf')
        close all;  tform1 = getRegistrationTransform(fixedMsi, unfixedMsi, 'surf');
        tform = tform1; 
    else 
        close all;  tform2 = getRegistrationTransform(fixedMsi, unfixedMsi, 'regconfig');
        tform = tform2;
    end
end 
newDims = [size(fixedMsi, 2), size(fixedMsi, 3)];
[recovered, unfixedWhiteReference, unfixedMask] = registerAllRelated(unfixedMsi, unfixedWhiteReference, unfixedMask, tform, newDims);

preprocUnfixed = recovered;
preprocFixed = fixedMsi; 

%% Identify crop areas 
[unfixedIaug, unfixedBbox, hasRoi] = getCropAreas(unfixedRoiCorners, roiNames, unfixedMask, unfixedWhiteReference);
[fixedIaug, fixedBbox, ~] = getCropAreas(fixedRoiCorners, roiNames, fixedMask, fixedWhiteReference);
    
metricsMel = zeros(length(mapMethods), length(metricsNames));
metricsHb = zeros(length(mapMethods), length(metricsNames));
roiMetricsMel = zeros(length(mapMethods), length(roiNames), length(metricsNames));
roiMetricsHb = zeros(length(mapMethods), length(roiNames), length(metricsNames));
imageNames = combineNameLists(mapMethods, tissueStates);

%% Chromophore Absorbance Maps
resultFixed = prepareMaps(true, preprocFixed, fixedMask, ID(fixedId), unfixedIaug, unfixedBbox, hasRoi);
resultUnfixed = prepareMaps(false, preprocUnfixed, unfixedMask, ID(unfixedId), fixedIaug, fixedBbox, hasRoi);

%% Plot ROIs
figTitle = strjoin({tissueStates{1}, 'vs', tissueStates{2}, strcat('(', lesion, ')')}, ' ' );
setSetting('plotName', fullfile(getSetting('savedir'), getSetting('common'), strcat(lesion, '_', 'rois.png')));
plotFunWrapper(1, @plotMontage, resultUnfixed.Iaug, resultFixed.Iaug, figTitle);

%% Plot Maps Collectively
plotFunWrapper(2, @plotMontageScaled, makeImageList(resultUnfixed.MelMaps, resultFixed.MelMaps), ...
    {unfixedMask, fixedMask}, imageNames, strcat(lesion, '_Full_Mel'));
plotFunWrapper(3, @plotMontageScaled, makeImageList(resultUnfixed.HbMaps, resultFixed.HbMaps), ...
    {unfixedMask, fixedMask}, imageNames, strcat(lesion, '_Full_Hb'));

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
    [metricsMel(i, :)] = getImageSimilarity(resultUnfixed.MelMaps{i}, resultFixed.MelMaps{i}, strcat(mapMethods{i}, '_mel'), lesion, tissueStates);
    [metricsHb(i, :)] = getImageSimilarity(resultUnfixed.HbMaps{i}, resultFixed.HbMaps{i}, strcat(mapMethods{i}, '_hb'), lesion, tissueStates);

    for k = 1:length(roiNames)
        if hasRoi(k)
            roiMetricsMel(i, k, :) = getImageSimilarity(resultUnfixed.RoiMelMaps{i, k}, resultFixed.RoiMelMaps{i, k}, ...
                strcat(mapMethods{i}, '_', roiNames{k}, '_mel'), lesion, tissueStates);
            roiMetricsHb(i, k, :) = getImageSimilarity(resultUnfixed.RoiHbMaps{i, k}, resultFixed.RoiHbMaps{i, k}, ...
                strcat(mapMethods{i}, '_', roiNames{k}, '_hb'), lesion, tissueStates);
        end
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

%% Save run
save(strrep(strcat(date, '_', lesion, '.mat'), ' ', '_')); 

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

function [roiCorners, lesion, regist] = getRoiCorners(idd, roiVals, txtVals)
idx = find(roiVals(:, 1) == idd);
roiCorners = {roiVals(idx, 2:5), roiVals(idx, 6:9), roiVals(idx, 10:13)};
lesion = txtVals{idx + 1, 1};
regist = strrep(txtVals{idx + 1, 15}, ' ', '');
end 

function [Iaug, bbox, hasRoi] = getCropAreas(roiCorners, roiNames, mask, white)
    roiColors = {'m', 'g', 'c'};
    roiMask = cell(numel(roiNames), 1);
    bbox = cell(numel(roiNames), 1);
    hasRoi = zeros(numel(roiNames));
    
    Iaug = white;
    for i = 1:numel(roiNames)
        [roiMask{i}, bbox{i}] = getBoundingBoxMask(roiCorners{i}, mask);
        hasRoi(i) = sum(isnan(bbox{i})) < 1;
        if hasRoi(i)
            Iaug = insertShape(Iaug, 'rectangle', bbox{i}, 'LineWidth', 5, 'Color', roiColors{i});
        end
    end
end 