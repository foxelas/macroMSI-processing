
%% Required settings

normType = 'divMacbeth';
removeBg = 'true';
fixedId = 16;
unfixedId = 17;
% msiType = 'extended';
msiType = 'adjusted';

%% Settings for experiment results 

roiNames = {'Hb', 'Norm', 'Mel'};
roiColors = {'m', 'g', 'c'};
mapMethods = {'Vasefi', 'Diebele', 'Kapsokalyvas'};
%{ 'vasefi', 'ding', 'ours', 'diebele', 'kapsokalyvas', 'kuzmina'};
roiCornerVals = delimread(fullfile(getSetting('systemdir'), getSetting('roiCornerFileName')), ',', 'num');
roiCornerVals = roiCornerVals.num;
idx = find(roiCornerVals(:, 1) == unfixedId);
roiCorners = {roiCornerVals(idx, 3:6), roiCornerVals(idx, 7:10), roiCornerVals(idx, 11:14)};
%roiCorners = {[316, 382, 242, 295], [159, 252, 167, 214], [398, 440, 83, 137]};

metricsNames = {'SSIM', 'NCC', 'HI', 'KLD', 'EMD', 'SSD', 'CC'};

%% Handle Fixed Image
[fixedMsi, fixedWhiteReference, fixedMask] = readAndNormalize(fixedId, ID, msiType, removeBg, normType);

%% Handle Unfixed Image
[unfixedMsi, unfixedWhiteReference, unfixedMask] = readAndNormalize(unfixedId, ID, msiType, removeBg, normType);

%% Register
setSetting('saveImages', true);
% close all;  tform1 = getRegistrationTransform(fixedMsi, unfixedMsi, 'surf');
close all;  tform2 = getRegistrationTransform(fixedMsi, unfixedMsi, 'regconfig');
newDims = [size(fixedMsi, 2), size(fixedMsi, 3)];
[recovered, unfixedWhiteReference, unfixedMask] = registerAllRelated(unfixedMsi, unfixedWhiteReference, unfixedMask, tform2, newDims);

preprocUnfixed = recovered;
preprocFixed = fixedMsi;


