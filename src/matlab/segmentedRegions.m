function [mask, maskAgreement, segmentMaskI] = segmentedRegions( files, coordinates, options, accTheta, regionRadius, thresVal, tfMean)
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
    if (nargin < 7) || isempty(tfMean)
        tfMean = false;
    end
         

    % Retrieve whole MSI
    [MSI, whiteReference, ~] = readMSI(files); %     I = squeeze(MSI(1,:,:,:));    
    g = permute(raw2msi(MSI, 'extended'), [2, 3, 1]);
    
    segmentMaskI = cell(ROIs,1);
    [m, n, bands] = size(g);
    bandWeight = 1 / bands;       
    totalMaskBinary = false(m,n);
    for roi = 1:ROIs
        maskAgreement = zeros(m, n); % RGB colored mask to show which MSI band resulted in which region
        mask = zeros(m, n); % Final mask for the region
        x = coordinates(roi, 1);
        y = coordinates(roi, 2);
        for i = 1:bands  
            I2 = squeeze(g(:,:,i));
            [~, maskTmp]= regionGrowing(I2, [y, x], thresVal, regionRadius, tfMean);
            maskAgreement = maskAgreement + bandWeight * maskTmp; % cat(3, bandColors(i,1) * maskTmp, bandColors(i,2) * maskTmp, bandColors(i,3) * maskTmp);
            mask = mask + bandWeight * maskTmp; 
        end

        mask = mask > accTheta; 
        if sum(mask(:)) < 9 %if the region is too small, add the neigborhood of the region seed pixel
            patchX = (x-2):(x+2);
            patchY = (y-2):(y+2);
            mask(patchY, patchX) = true;
        end
        segmentMaskI{roi} = mask;   
        totalMaskBinary = totalMaskBinary + mask;
    end
    
    if (options.showImages)
        plots('segmentation', 1, 'Image', whiteReference, 'Overlay', maskAgreement, ...
            'AdditionalImage',  whiteReference + mask, 'Coordinates', coordinates, 'SaveOptions', options.saveOptions);
    end
end
