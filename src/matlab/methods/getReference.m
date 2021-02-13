function [reference] = getReference(systemdir, h, w)
%     GETREFERENCE returns white reference image
%
%     Usage:
%     reference = getReference(systemdir, h, w)


infile = fullfile(systemdir, 'infiles', getSetting('whiteReference'));
load(infile, 'raw');
[~, height, width, ~] = size(raw);

if (nargin < 2)
    h = height;
    w = width;
end

halvingLowFun = @(x) floor(x/2);
xCenter = halvingLowFun(height);
yCenter = halvingLowFun(width);
hStart = floor(xCenter-h/2);
wStart = floor(yCenter-w/2);

reference = raw(:, hStart:h+hStart-1, wStart:w+wStart-1, :);
setSetting('saveImages', false);
plots(1, @plotMSI, reference, false);
setSetting('saveImages');
end
