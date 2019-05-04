function [segmentMSI, segmentWhite, segmentDark, segmentMask, segmentMaskI] = segmentMSIRegion( files, coordinates, options, accTheta, regionRadius, thresVal)
%%SEGMENTMSIREGION applies region growing on every channel of the MSI for
%%seed position coordinates = [x,y], based on region agreement threshold 'accTheta' and 
%%radius 'regionRadius'
%The resulting mask is decided per agreement of accTheta * 1/NumberOfChannels

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
    [MSI, whiteReference, darkReference] = readMSI(files); %     I = squeeze(MSI(1,:,:,:));    
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

        if (options.showImages)
            currentOptions = options.saveOptions;
            currentOptions.plotName = options.saveOptions.plotName{roi};
            plots('segmentation', 1, 'Image', whiteReference, 'Overlay', maskAgreement, ...
                'AdditionalImage',  whiteReference + mask, 'Coordinates', [x,y], 'SaveOptions', currentOptions);
        end

        if ~isempty(whiteReference)
            segmentWhite{roi} = whiteReference(patchY, patchX,:);
        end
        if ~isempty(darkReference)
            segmentDark{roi} = darkReference(patchY, patchX,:);
        end
        
    end
end
