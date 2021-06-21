function normMsi = getNormalizedMsi(raw, normType, msiType)
%     GETNORMALIZEDMSI returns the normalized msi

%     Input arguments
%     raw: the RAW multispectral image
%     normType: the normalization type ['none', 'divByMax', 'divMacbeth']
%     msiType: the msi construction type  {'green', 'rms', 'adjusted',
%     'extended', 'unchanged', 'max'}
%
%     Output arguments
%     normMsi: the normalized msi
%
%     Usage:
%     normMsi = getNormalizedMsi(raw, 'divByMax');
%     normMsi = getNormalizedMsi(raw, 'divByMax', 'extended');

if nargin < 3
    msiType = 'extended';
end

msi = raw2msi(raw, msiType);

switch normType
    case 'none'
        normMsi = msi;
    case 'divByMax'
        normRaw = raw ./ max(raw(:));
        normMsi = raw2msi(normRaw, msiType);
    case 'subAvg'
        normMsi = getDifMSI(msi, 'toAverage');
    case 'divMacbeth'
        [~, height, width] = size(msi);
        reference = getReferenceFromMacbeth(height, width);
        reference = raw2msi(reference, msiType);
        normMsi = msi ./ reference;
end
end
