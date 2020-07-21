function [result] = getReferenceFromMacbeth(h, w)
%     GETREFERENCEFROMMACBETH returns white reference image 
% 
%     Usage: 
%     reference = getReferenceFromMacbeth(h, w)
load(getSetting('whiteReferenceMacbeth'), 'referenceWhite');

[ channels, ~, rgbDim] = size(referenceWhite);
result = reshape(repmat(referenceWhite, [1, h * w, 1]), [channels, h, w, rgbDim]);

% oldSetting = getSetting('saveImages');
% setSetting('saveImages', false);
% plotFunWrapper(1, @plotMSI, result);
% setSetting('saveImages', oldSetting);

end 