function [segmentMSI, segmentWhite, segmentDark, segmentMask, segmentMaskI] = readMSI(files, coordinates, width, height, options, fc)
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
    
    load('saved parameters\color_correction.mat', 'illuminant_gw1');

    extraImages = 0;
    [hasWhiteReference, idx] = ismember(1, fc);
    whiteReference = [];
    if hasWhiteReference
        whiteReference = chromadapt(imread(files{idx}), illuminant_gw1, 'ColorSpace', 'linear-rgb');
        whiteReference = im2double(whiteReference);
        [imHeight, imWidth, ~] = size( whiteReference );
        extraImages = extraImages + 1;
    end

    [hasDarkReference, idx] = ismember(0, fc);
    darkReference = [];
    if hasDarkReference
        darkReference = chromadapt(imread(files{idx}), illuminant_gw1, 'ColorSpace', 'linear-rgb');
        darkReference = im2double(darkReference);
        extraImages = extraImages + 1;
    end

    MSIbands = length(files) - extraImages;
    MSI = zeros(MSIbands, imHeight, imWidth, 3);
    for k = 1:MSIbands
        msi = chromadapt(imread(files{k+extraImages}), illuminant_gw1, 'ColorSpace', 'linear-rgb');
        MSI(k, :, :, :) = im2double(msi);
    end
    
    if (modeAll)
        segmentWhite = whiteReference;
        segmentDark = darkReference;
        segmentMSI = MSI;
        segmentMaskI = ones(imHeight, imWidth);
        segmentMask = ones(imHeight, imWidth);
        
    else
        segmentWhite = cell(ROIs,1);
        segmentDark = cell(ROIs,1);
        segmentMSI = cell(ROIs,1);
        segmentMask = cell(ROIs,1);
        segmentMaskI = cell(ROIs,1);
        
        for roi = 1:ROIs
            
            x = coordinates(roi, 1);
            y = coordinates(roi, 2);
            [patchMask, maskI] = makeMasks(imHeight, imWidth, x, y, width, height);
            
            if hasWhiteReference
                segmentWhite{roi} = whiteReference(y:(y + height - 1), x:(x + width - 1), :);
                if isfield(options, 'showImages') && (options.showImages)
                    currentOptions = options.saveOptions;
                    currentOptions.plotName = options.saveOptions.plotName{roi};
                    plots('cropped', 1, 'Image', whiteReference + maskI, 'Coordinates', [x,y], 'SaveOptions', currentOptions);
                end
            end
            
            if hasDarkReference
                segmentDark{roi} = darkReference(y:(y + height - 1), x:(x + width - 1), :);
            end
            
            segmentMSI{roi} = MSI(:, y:(y + height - 1), x:(x + width - 1), :);
            segmentMask{roi} = patchMask;
            segmentMaskI{roi} = maskI;

        end  
        
        if (ROIs == 1) %if it it reads only one ROI, there is no need for cell array
            segmentMSI = cell2mat(segmentMSI);
            segmentMask = cell2mat(segmentMask);
            segmentMaskI = cell2mat(segmentMaskI);
            segmentDark = cell2mat(segmentDark);
            segmentWhite = cell2mat(segmentWhite);        
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
