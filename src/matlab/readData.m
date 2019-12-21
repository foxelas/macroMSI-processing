%% LOAD DATA & INITIALIZE

disp('Initialization.')
load(fullfile(options.systemdir, 'system.mat')); % camera system parameters
load(fullfile(options.systemdir, 'data.mat')); % image data
load(fullfile(options.systemdir, 'ID.mat')); % image data id and info struct
%load(fullfile(options.systemdir, 'precomputedParams.mat')); %pre-set parameters

if (contains(dataset, 'bright'))
    ID = IDbrights;
else 
    ID = IDdarks;
end

wavelengthN = size(sensitivity, 1);
wavelength = linspace(380, 780, wavelengthN);
originalChannels = 7;
outputFolderMap = getOutputDirectoryMap();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msiN = length(ID);
if (options.tryReadData)
    %w = warning('off', 'all');   
    if  ~isfile(generateName('matfilein', options))
        dateCreated = datetime();
        save(generateName('matfilein', options), 'dateCreated');
    end

    %% Read raw spectra 
    [Spectra,CompleteSpectra,SpectraNames] = bulkReadSpectra(ID, data, wavelengthN, options);
    
    %% Read MSI
    fprintf('Reading MSI data according to ID file.\n');

    groups = findgroups([ID.MsiID]);
    groupN = max(groups);
     
    msiN = length(ID);  
    g = 0;
    
    for k = 1:msiN 
        idd = ID(k);
        if g ~= groups(k)
            g = groups(k);

            gIdxs = find(groups == g);
            gMembers = ID(gIdxs);
            coordinates = [gMembers.Originx; gMembers.Originy];
            files = {data([data.MsiID] == idd.MsiID).File};
            [raw, whiteReference,darkReference] = readMSI(files);
            specimenMask = removeBackground(whiteReference, idd, options);
            
            infile = mkdir_custom(fullfile(options.systemdir, 'infiles', strcat('group_', num2str(idd.Group), '.mat')));
            save(infile, 'raw', 'whiteReference', 'darkReference', 'specimenMask');
        end
        
        segmentMask = segment(raw, idd, whiteReference,...
            specimenMask, options, []);
        segmentMasksForFeatureExtraction = segment(raw, idd, whiteReference, specimenMask, ...
            options, 50);
        poiName = strjoin({num2str(idd.Index), idd.Sample, idd.Type, idd.T, idd.Label}, '_');
        
        [r, c] = find(segmentMasksForFeatureExtraction);
        patchY = min(r):max(r);
        patchX = min(c):max(c);
        poiRAW = raw(:, patchY, patchX, :);
        poiWhite = whiteReference(patchY, patchX, :);
        poiSegmentMask = segmentMasksForFeatureExtraction(patchY, patchX);
        roiSeeds = [ID(k).Originx - min(c) + 1, ID(k).Originy - min(r) + 1];
        measuredSpectrum = Spectra(k,:);

        infile = mkdir_custom(fullfile(options.systemdir, 'infiles', strcat('poi_', num2str(idd.Index), '.mat')));
        save(infile, 'segmentMask', 'segmentMasksForFeatureExtraction', 'poiName',...
            'poiRAW', 'poiSegmentMask', 'roiSeeds', 'measuredSpectrum', 'poiWhite');

    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Precompute correlation smoothing Matrices 
    computeSmoothingMatrices(Spectra, SpectraNames, wavelengthN, ID, options);
    load(precomputedFile);

    %% IlluminationXSensitivity
    Hgreen = illumXsensitivity(illumination, sensitivity, 'green');
    Hadjusted = illumXsensitivity(illumination, sensitivity, 'adjusted');
    Hrms = illumXsensitivity(illumination, sensitivity, 'rms');
    Hextended = illumXsensitivity(illumination, sensitivity, 'extended');

    save(precomputedFile, 'Hgreen', 'Hadjusted', 'Hrms', 'Hextended', '-append');

    %% Create coefficient table
    createCoefficientTable(ID, CompleteSpectra, options);
    load( precomputedFile, 'Coefficients');
    
    %% Search for reference coefficients 
    ID = selectCoeffIndex(options); 
end
    
fprintf('Finished initialization.\n\n')
fprintf('Let''s get down to business :muscle:\n\n');
% LOAD DATA & INITIALIZE ends   
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H] = illumXsensitivity(E, S, pixelValueSelectionMethod)
%illumXsensitivity Computes the system matrix H using illumination and
%sensitivity
%   Input args
%   E: illumination , spectrumdim x MSdim e.g. 401x7
%   S: camera sensitivity, spectrumdim x filterdim, e.g. 401x3, RGB
%   pixelValueSelectionMethod: either of {'green', 'rms', 'adjusted', 'extended', 'unchanged'}
%   Output args
%   H: system matrix, MSdim x spectrumdim e.g. 7x401

    if (size(S, 2) > 3)
        if strcmp(pixelValueSelectionMethod, 'extended')
            pixelValueSelectionMethod = 'fullbandext';
        else
            pixelValueSelectionMethod = 'fullband';
        end
    end
    
    if strcmp(pixelValueSelectionMethod, 'rgb')
        E = E(:, 8); % white light
    else
        E = E(:, 1:7); % narrow bandwidth lights
    end

    switch pixelValueSelectionMethod
        case 'green'
            H = E' * diag(S(:, 2));

        case 'rms'
            H = E' * diag(sqrt(sum(S.^2, 2)));

        case 'adjusted'
            H(1:2, :) = E(:, 1:2)' * diag(S(:, 3)); %  blue range
            H(3:5, :) = E(:, 3:5)' * diag(S(:, 2)); %  green range
            H(6:7, :) = E(:, 6:7)' * diag(S(:, 1)); %  red range

        case 'extended'
            H(1:3, :) = E(:, 1:3)' * diag(S(:, 3)); % extended blue range
            H(4:7, :) = E(:, 3:6)' * diag(S(:, 2)); % extended green range
            H(8:9, :) = E(:, 6:7)' * diag(S(:, 1)); % extended red range

        case 'unchanged' %estimation using RGB images
            H = zeros(3, length(S));
            lightBands = [6, 7; 3, 5; 1, 2]; %RGB light bands
            for k = 1:3
                Eadj = rms(E(:, lightBands(k, 1):lightBands(k, 2)), 2)'; % average the light at certain wavelength bands
                H(k, 1:length(S)) = Eadj .* S(:, k)';
            end

        case 'rgb'
            H(1, :) = E' * diag(S(:, 3)); % extended blue range
            H(2, :) = E' * diag(S(:, 2)); % extended green range
            H(3, :) = E' * diag(S(:, 1)); % extended red range

        case 'fullband'
            H = (E .* S)';

        case 'fullbandext'
            E = [E(:, 1:3), E(:, 3:6), E(:, 6:7)];
            S = [S(:, 1:3), S(:, 3:6), S(:, 6:7)];
            H = (E .* S)';

        otherwise
            error('Unexpected pixelValueSelectionMethod. Could not compute system matrix.')
    end

end

