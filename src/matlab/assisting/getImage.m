function [msi, whiteReference, specimenMask, height, width, channels] = getImage(k, options, msiType, removebg)

    if nargin < 3 
        msiType = 'max'; %'extended'; % 'max';
    end 
    
    if nargin < 4 
        removebg = true; 
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

    options.saveOptions.saveImages = false; 
    plotMSI(msi, 1, options.saveOptions);
end
