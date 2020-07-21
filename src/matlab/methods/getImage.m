function [msi, whiteReference, specimenMask, height, width, channels] = getImage(k, msiType, removebg, isColumnImage, normType, tform, newDims)
%     GETIMAGE returns the msi and other useful parameters
% 
%     Input arguments
%     k: the index of the requested image file f
%     msiType: the msi construction type  {'green', 'rms', 'adjusted',
%     'extended', 'unchanged', 'max'}
%     removeBg: the boolean variable showing whether background pixels should
%     be set to 0 (if TRUE) or background pixels are left as is (if FALSE)
%     isColumnImage: the boolean variable showing whether the each MSI band
%     subimage should be covnerted to a column vector (if TRUE) or remain as
%     is (if FALSE)
%     normType: the normalization type ['none', 'divByMax', 'divMacbeth']
%     tform: translation transform for rotation and scaling 
%     newDims: the new dimensions for the translated image
% 
%     Output arguments
%     msi: the multispectral image
%     whiteReference: the respective RGB image
%     specimenMask: the specimen binary mask
%     height: the 2D channel subimage height
%     width: the 2D channel subimage width
%     channels: the number of channels
%
%     Usage: 
%     [msi, whiteReference, specimenMask, height, width, channels] = getImage(k)
%     [msi, whiteReference, specimenMask, height, width, channels] = getImage(k, msiType)
%     [msi, whiteReference, specimenMask, height, width, channels] = getImage(k, msiType, removebg)
%     [msi, whiteReference, specimenMask, height, width, channels] = getImage(k, msiType, removebg, isColumnImage)

if nargin < 2
    msiType = 'max'; %'extended'; % 'max';
end

if nargin < 3
    removebg = true;
end

if nargin < 4
    isColumnImage = false;
end

if nargin < 5 
    normType  = 'none';
end 

if nargin < 6 
    tform = [];
    newDims = [];
end 

infile = fullfile(getSetting('systemdir'), 'infiles', strcat('group_', num2str(k), '.mat'));
load(infile, 'raw', 'whiteReference', 'specimenMask');

if ~isempty(tform) && ~isempty(newDims) 
    raw = registerImage(raw, tform, newDims);
    whiteReference = permute(registerImage(permute(whiteReference, [3, 1, 2]), tform, newDims), [2, 3, 1]);
    specimenMask = registerImage(specimenMask, tform, newDims);
end 

[~, height, width, ~] = size(raw);
msi = getNormalizedMsi(raw, normType, msiType);
[channels, ~, ~] = size(msi);


if removebg
    foregroundMask = permute(repmat(double(specimenMask), 1, 1, channels), [3, 1, 2]);
    msi = bsxfun(@times, msi, foregroundMask);
end

if isColumnImage
    columns = reshape(msi, channels, width*height)'; %30k pixels x 7 variable
    columnsWhite = reshape(permute(whiteReference, [3, 1, 2]), 3, width*height)'; %30k pixels x 3 variable
    fgColumn = reshape(specimenMask, 1, width*height);

    if removebg
        columns = columns(fgColumn, :);
        columnsWhite = columnsWhite(fgColumn, :);
    end
    msi = columns;
    whiteReference = columnsWhite;
    specimenMask = fgColumn;
else
    setSetting('saveImages', false);
    plotFunWrapper(1, @plotMSI, msi, false);
    getSetting('saveImages');
end
end
