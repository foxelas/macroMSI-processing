function [segments] = readMSI(files, coordinates, width, height, options, fc)
%Reads the MSI image file and contains it in a matrix g
% files: the subimages related to the MSI
% x: origin x to read image box
% y: origin y to read image box
% width: width of image box to be read (x axis)
% height: height of image box to be read (y axis)

    if ~exist('coordinates', 'var') 
        modeAll = true; 
        
    elseif isempty(coordinates) 
        modeAll = true;   
        
    else
        modeAll = false;
        if (size(coordinates, 1) == 2)
            coordinates = coordinates';
        end
        [ROIs, coor] = size(coordinates); % m is the number of ROIs to be read from the MSI, n is 2 for (x,y)
        if (coor ~= 2)
            error('Coordinates should consist of 2 dimensions.');
        end
    end

    if (nargin < 3)
        width = 5;
    end

    if (nargin < 4)
        height = 5;
    end

    % modeAll: if true, reads the whole image
    if (width == 0) || (height == 0)
        error('No region to be read.')
    end

    if (nargin < 5) || isempty(options)
        options.saveOptions.plotName = '../output/cropped/';
    end

    if (nargin < 6) || isempty(fc)
        fc = [1, 450, 465, 505, 525, 575, 605, 630];
    end   
    
    segments = struct('msi', [], 'whiteI', [], 'darkI', [], 'patchMask', [], 'maskI', []);

    extraImages = 0;
    [hasWhiteReference, idx] = ismember(1, fc);
    whiteReference = [];
    if hasWhiteReference
        whiteReference = imread(files{idx});
        [imHeight, imWidth, ~] = size( whiteReference );
        extraImages = extraImages + 1;
    end

    [hasDarkReference, idx] = ismember(0, fc);
    darkReference = [];
    if hasDarkReference
        darkReference = imread(files{idx});
        extraImages = extraImages + 1;
    end

    MSIbands = length(files) - extraImages;
    MSI = zeros(MSIbands, imHeight, imWidth, 3);
    for k = 1:MSIbands
        MSI(k, :, :, :) = imread(files{k+extraImages});
    end
        
    if (modeAll)
        segments(1).whiteI = whiteReference;
        segments(1).darkI = darkReference;
        segments(1).msi = MSI;
        segments(1).maskI = ones(imHeight, imWidth);
        
    else
    
        for roi = 1:ROIs
            
            x = coordinates(roi, 1);
            y = coordinates(roi, 2);
            
            if hasWhiteReference
                segments(roi).whiteI = im2double( whiteReference(y:(y + height - 1), x:(x + width - 1), :) );
                showCroppedSection( whiteReference, '', x, y, strrep([files{1} ,' ', num2str(id)], '_', ' '),  [options.saveOptions.plotName, '_cropped.jpg']);
            end
            
            if hasDarkReference
                segments(roi).darkI = im2double( darkReference(y:(y + height - 1), x:(x + width - 1), :) );
            end
            
            segments(roi).msi = im2double( MSI(:, y:(y + height - 1), x:(x + width - 1), :) );
            [segments(roi).patchMask, segments(roi).maskI] = makeMasks( imHeight, imWidth, x, y, width, height);

        end  
        
    end

end

function  [patchMask, maskI] = makeMasks(imHeight, imWidth, x, y, width, height)
    
    patchMask = ones(height, width);
    maskI = zeros(imHeight, imWidth);
    if (x > 0) && (y > 0)
        maskI(y:(y + height - 1), x:(x + width - 1)) = 1;
    end
    
end
