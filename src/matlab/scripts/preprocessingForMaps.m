
%% Required settings
% msiType = 'extended';
% normType = 'divMacbeth';
% removeBg = 'true';
% fixedId = 16;
% unfixedId = 17;

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