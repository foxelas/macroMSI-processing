
%% Required settings

normType = 'divMacbeth';
removeBg = 'true';
fixedId = 16;
unfixedId = 17;
% msiType = 'extended';
msiType = 'adjusted';

%% Settings for experiment results 

%{ 'vasefi', 'ding', 'ours', 'diebele', 'kapsokalyvas', 'kuzmina'};
roiCornerVals = delimread(fullfile(getSetting('systemdir'), getSetting('roiCornerFileName')), ',', {'num', 'text'});
diseaseVals = roiCornerVals.text; 
roiCornerVals = roiCornerVals.num;
[unfixedRoiCorners, ~] = getRoiCorners(unfixedId, roiCornerVals, diseaseVals); 
[fixedRoiCorners, lesion] = getRoiCorners(fixedId, roiCornerVals, diseaseVals); 
%roiCorners = {[316, 382, 242, 295], [159, 252, 167, 214], [398, 440, 83, 137]};

tissueStates = {'unfixed', 'fixed'}; 
metricsNames = {'SSIM', 'NCC', 'HI', 'KLD', 'EMD', 'SSD', 'CC'};

metricsMel = zeros(length(mapMethods), length(metricsNames));
metricsHb = zeros(length(mapMethods), length(metricsNames));
roiMetricsMel = zeros(length(mapMethods), length(roiNames), length(metricsNames));
roiMetricsHb = zeros(length(mapMethods), length(roiNames), length(metricsNames));
imageNames = combineNameLists(mapMethods, tissueStates);


% %% Handle Fixed Image
% [fixedMsi, fixedWhiteReference, fixedMask] = readAndNormalize(fixedId, ID, msiType, removeBg, normType);
% 
% %% Handle Unfixed Image
% [unfixedMsi, unfixedWhiteReference, unfixedMask] = readAndNormalize(unfixedId, ID, msiType, removeBg, normType);
% 
% %% Register
% setSetting('saveImages', true);
% % close all;  tform1 = getRegistrationTransform(fixedMsi, unfixedMsi, 'surf');
% close all;  tform2 = getRegistrationTransform(fixedMsi, unfixedMsi, 'regconfig');
% newDims = [size(fixedMsi, 2), size(fixedMsi, 3)];
% [recovered, unfixedWhiteReference, unfixedMask] = registerAllRelated(unfixedMsi, unfixedWhiteReference, unfixedMask, tform2, newDims);
% 
% preprocUnfixed = recovered;
% preprocFixed = fixedMsi;

function [roiCorners, lesion] = getRoiCorners(idd, roiVals, disVals)
idx = find(roiVals(:, 1) == idd);
roiCorners = {roiVals(idx, 2:5), roiVals(idx, 6:9), roiVals(idx, 10:13)};
lesion = disVals{idx + 1};
end 