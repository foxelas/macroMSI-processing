function [msi, whiteReference, specimenMask, height, width, channels] = getImage(k, options, msiType, removebg, isColumnImage)
%% FUNCTION getImage returns the msi and other useful parameters 
%
%   Input arguments 
%   k: the index of the requested image file 
%   options: the running options 
%   msiType: the msi construction type  {'green', 'rms', 'adjusted', 
%   'extended', 'unchanged', 'max'}
%   removeBg: the boolean variable showing whether background pixels should 
%   be set to 0 (if TRUE) or background pixels are left as is (if FALSE)
%   isColumnImage: the boolean variable showing whether the each MSI band
%   subimage should be covnerted to a column vector (if TRUE) or remain as
%   is (if FALSE)
%   
%   Output arguments 
%   msi: the multispectral image 
%   whiteReference: the respective RGB image 
%   specimenMask: the specimen binary mask
%   height: the 2D channel subimage height 
%   width: the 2D channel subimage width 
%   channels: the number of channels
%

    if nargin < 3 
        msiType = 'max'; %'extended'; % 'max';
    end 
    
    if nargin < 4 
        removebg = true; 
    end 
    
    if nargin  < 5 
        isColumnImage = false;
    end 
    
    infile = fullfile(options.systemdir, 'infiles', strcat('group_', num2str(k), '.mat'));
    load(infile, 'raw', 'whiteReference', 'specimenMask');
    
    [~, height, width, ~] = size(raw);
    msi = raw2msi(raw, msiType);
    [channels, ~, ~] = size(msi);

    if removebg 
        foregroundMask = permute(repmat(double(specimenMask), 1, 1,  channels), [3 1 2]);   
        msi = bsxfun(@times, msi, foregroundMask);
    end 
    
    if isColumnImage
        columns = reshape(msi, channels, width * height)'; %30k pixels x 7 variable   
        columnsWhite = reshape(permute(whiteReference, [3 1 2]), 3, width * height)'; %30k pixels x 3 variable     
        fgColumn = reshape(specimenMask, 1, width*height);

        if removebg
            columns = columns(fgColumn, :);
            columnsWhite = columnsWhite(fgColumn, :);
        end
        msi = columns; 
        whiteReference = columnsWhite;
        specimenMask = fgColumn; 
    else 
        options.saveOptions.saveImages = false; 
        plotMSI(msi, 1, options.saveOptions);
    end
end
