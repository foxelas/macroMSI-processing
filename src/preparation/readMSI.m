function [segments] = readMSI(files, coordinates, width, height, options, fc)
%Reads the MSI image file and contains it in a matrix g
% files: the subimages related to the MSI
% coordinates: [x,y] origin to read image box
% width: width of image box to be read (x axis)
% height: height of image box to be read (y axis)
% options: running configurations
% fc: the frequencies of the MSI 

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
    
    segments = struct('MSI', [], 'whiteReference', [], 'darkReference', [], 'patchMask', [], 'maskI', []);

    extraImages = 0;
    [hasWhiteReference, idx] = ismember(1, fc);
    whiteReference = [];
    if hasWhiteReference
        whiteReference = im2double(imread(files{idx}));
        [imHeight, imWidth, ~] = size( whiteReference );
        extraImages = extraImages + 1;
    end

    [hasDarkReference, idx] = ismember(0, fc);
    darkReference = [];
    if hasDarkReference
        darkReference = im2double(imread(files{idx}));
        extraImages = extraImages + 1;
    end

    MSIbands = length(files) - extraImages;
    MSI = zeros(MSIbands, imHeight, imWidth, 3);
    for k = 1:MSIbands
        MSI(k, :, :, :) = im2double(imread(files{k+extraImages}));
    end
    
    if (modeAll)
        segments(1).whiteReference = whiteReference;
        segments(1).darkReference = darkReference;
        segments(1).MSI = MSI;
        segments(1).maskI = ones(imHeight, imWidth);
        
    else
        
        for roi = 1:ROIs
            
            x = coordinates(roi, 1);
            y = coordinates(roi, 2);
            [patchMask, maskI] = makeMasks(imHeight, imWidth, x, y, width, height);
            
            if hasWhiteReference
                segments(roi).whiteReference = whiteReference(y:(y + height - 1), x:(x + width - 1), :);
                if (options.showImages)
                    currentOptions = options.saveOptions;
                    currentOptions.plotName = options.saveOptions.plotName{roi};
                    plots('cropped', 1, 'Image', whiteReference + maskI, 'Markers', 'r*', 'Coordinates', [x,y], 'SaveOptions', currentOptions);
                end
            end
            
            if hasDarkReference
                segments(roi).darkReference = darkReference(y:(y + height - 1), x:(x + width - 1), :);
            end
            
            segments(roi).MSI = MSI(:, y:(y + height - 1), x:(x + width - 1), :);
            segments(roi).patchMask = patchMask;
            segments(roi).maskI = maskI;

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
