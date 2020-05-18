function [segmentMSI, segmentWhite, segmentDark, segmentMask, segmentMaskI] = segmentMSIRegion( files, coordinates, accTheta, regionRadius, thresVal, plotNames)
%%segmentMSIRegion Read an MSI segment from raw RGB subimages
%
% [segmentMSI, segmentWhite, segmentDark, segmentMask, segmentMaskI] = 
%     segmentMSIRegion( files, coordinates, accTheta, 
%     regionRadius, thresVal)
% Reads subimages with given filenames and loads the area obtained after 
% region growing from original point at coordinates to a 4D matrix. 
% Additionally, produces the respective white and dark RGB images and 
% binary masks of the read area.
% 
% Inputs:
% files - filenames of subimages related to the MSI
% coordinates - [x,y] of the seed pixel for region growing 
% accTheta - agreement threshold for the MSI channel 
% regionRadius - maximum radius of the grown region
% thresVal - pixel intensity threshold 
% 
% Outputs:
% segmentMSI - the 4D raw MSI created from RGB subimages
% segmentWhite - the respective RGB image under white light
% segmentDark - the respective RGB image under no illumination
% segmentMask - binary mask of the points of interest inside the bounding box
% segmentMaskI - binary masks of the points of interest inside the entire 
%     dimensions of the MSI
%     
% Usage:
% [segmentMSI, segmentWhite, segmentDark, segmentMask, segmentMaskI] = 
%     readMSI(files, coordinates, width, height,  fc)
% Reads the area of the image contained in bounding box with upper 
% left corner at [coordinates], with dimensions [width], [height]
% segmentMSI = readMSI(files)
% Reads the entire MSI image
%
    if (size(coordinates, 1) == 2)
        coordinates = coordinates';
    end
    [ROIs, coor] = size(coordinates); % m is the number of ROIs to be read from the MSI, n is 2 for (x,y)
    if (coor ~= 2)
        error('Coordinates should consist of 2 dimensions.');
    end
    
    if (nargin < 4) || isempty(accTheta)
        accTheta = 0.7;
    end
    if (nargin < 5) 
        regionRadius = 36;
    end

    if (nargin < 6) || isempty(thresVal)
        thresVal = [];
    end

    % Retrieve whole MSI
    [MSI, whiteReference, darkReference] = readMSI(files); 
    g = permute(raw2msi(MSI, 'extended'), [2, 3, 1]);
    
    [m, n, bands] = size(g);
    bandWeight = 1 / bands;       
    
    segmentWhite = cell(ROIs,1);
    segmentDark = cell(ROIs,1);
    segmentMSI = cell(ROIs,1);
    segmentMask = cell(ROIs,1);
    segmentMaskI = cell(ROIs,1);

    for roi = 1:ROIs
        
        maskAgreement = zeros(m, n);  % RGB colored mask to show which MSI band resulted in which region
        mask = zeros(m, n); % Final mask for the region 
        x = coordinates(roi, 1);
        y = coordinates(roi, 2);
        for i = 1:bands  
            I2 = squeeze(g(:,:,i));
            [~, maskTmp]= regionGrowing(I2, [y, x], thresVal, regionRadius);
            maskAgreement = maskAgreement + bandWeight * maskTmp;
            mask = mask + bandWeight * maskTmp; 
        end

        mask = mask > accTheta; 
        if sum(mask(:)) < 9 %if the region is too small, add the neigborhood of the region seed pixel
            patchX = (x-2):(x+2);
            patchY = (y-2):(y+2);
            mask(patchY, patchX) = 1;
        end
        segmentMaskI{roi} = mask;
        [r, c] = find(mask);
        patchY = min(r):max(r);
        patchX = min(c):max(c);
        segmentMSI{roi} = MSI(:, patchY, patchX, :);
        segmentMask{roi} = mask(patchY, patchX);

        if getSetting('showImages')
            setOption('plotName', plotNames{roi}); 
            plotSegmentation(whiteReference, maskAgreement, whiteReference + mask, [x,y], 1);
        end

        if ~isempty(whiteReference)
            segmentWhite{roi} = whiteReference(patchY, patchX,:);
        end
        if ~isempty(darkReference)
            segmentDark{roi} = darkReference(patchY, patchX,:);
        end
        
    end
end
