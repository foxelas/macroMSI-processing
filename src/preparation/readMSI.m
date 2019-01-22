function [MSI, whiteReference, darkReference] = readMSI(files, x, y, width, height, options, fc)
%Reads the MSI image file and contains it in a matrix g
% files: the subimages related to the MSI
% x: origin x to read image box
% y: origin y to read image box
% width: width of image box to be read (x axis)
% height: height of image box to be read (y axis)

if isempty(x) || isempty(y)
        modeAll = true;      
elseif (x < 0) || (y < 0)
    error('Negative origin points.')
else
        modeAll = false;
end

if (nargin < 4)
    width = 5;
end

if (nargin < 5)
    height = 5;
end

% modeAll: if true, reads the whole image
if (width == 0) || (height == 0)
    error('No region to be read.')
end

if (nargin < 6) || isempty(options)
    options.saveOptions.plotName = '../output/cropped/q';
end

if (nargin < 7) || isempty(fc)
    fc = [1, 450, 465, 505, 525, 575, 605, 630];
end

    function [SMSI] = readSubimage(loc)
        ImFull = imread(files{loc});
        if ~(modeAll)
            %             fprintf('now at %d %d\n', x, y)
            I = ImFull(y:(y + height - 1), x:(x + width - 1), :);
        else
            I = ImFull;
        end
        
        %to check that we read the correct region
        if (loc == 1 && fc(1) == 1 && ~modeAll) || (loc == 2 && fc(1) ~= 1 && ~modeAll)
            showCroppedSection( ImFull, '', x, y, strrep([files{loc} ,' ', num2str(id)], '_', ' '),  [options.saveOptions.plotName, '_cropped.jpg']);
        end
        
        SMSI = im2double(I);
    end

%  fig1 = figure(1);
extraImages = 0;
[hasWhiteReference, idx] = ismember(1, fc);
if hasWhiteReference
    whiteReference = readSubimage(idx);
    files{idx} = [];
    extraImages = extraImages + 1;
else
    whiteReference = [];
end

[hasDarkReference, idx] = ismember(0, fc);
if hasDarkReference
    darkReference = readSubimage(idx);
    files{idx} = [];
    extraImages = extraImages + 1;
else
    darkReference = [];
end

MSIbands = length(files) - extraImages;
for k = 1:MSIbands
    MSI(k, :, :, :) = readSubimage(k+extraImages);
end

end
