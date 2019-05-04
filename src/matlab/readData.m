%% LOAD DATA & INITIALIZE

disp('Initialization.')
load(fullfile(options.systemdir, 'system.mat')); % camera system parameters
load(fullfile(options.systemdir, 'data.mat')); % image data
load(fullfile(options.systemdir, 'ID.mat')); % image data id and info struct

wavelengthN = size(sensitivity, 1);
wavelength = linspace(380, 780, wavelengthN);
step = ceil(length(380:780) / wavelengthN);
wavelengthIdxs = 1:step:length(380:780);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msiN = length(ID);
if (options.tryReadData)
    w = warning('off', 'all');   
    fprintf('Reading spectral data according to ID file.\n');

    if  ~isfile(generateName(options, 'matfilein'))
        dateCreated = datetime();
        save(generateName(options, 'matfilein'), 'dateCreated');
    end

    %% Read raw spectra 
    [uniqueSpectraNames, uniqueSpectraIdxs, uniqueSpectraIdxsInID] = unique(strcat({ID.SpectrumFile}, {ID.T}));
    idd = ID(uniqueSpectraIdxs);
    specN = length(uniqueSpectraIdxs);
    completeUniqueSpectra = zeros(specN, 401);

    for i = 1:specN
        % read raw measured spectrum
        rawSpectrum = readSpectrum(idd(i).SpectrumFile, idd(i).T);
        % read raw white measured spectrum of the reference surface
        referenceSpectrum = readSpectrum(char(strcat(data(idd(i).RgbID).Sample, '\', 'white.csv')));

        if abs(rawSpectrum-referenceSpectrum) < 0.000001
            error('Measurement is same as white.')
        end
        completeUniqueSpectra(i,:) = rawSpectrum ./ referenceSpectrum;  
    end
    uniqueSpectra = completeUniqueSpectra(:,wavelengthIdxs);
    clear('idd');

    %% Save Spectra Struct
    Spectra = zeros(msiN, wavelengthN);
    CompleteSpectra = zeros(msiN, 401);
    SpectraNames = cell(msiN, 1);
    for i = 1:msiN            
        [~, csv] =  generateName(options, 'csv', ID(i));
        Spectra(i,:) = uniqueSpectra(uniqueSpectraIdxsInID(i),:);
        CompleteSpectra(i,:) = completeUniqueSpectra(uniqueSpectraIdxsInID(i),:);
        SpectraNames{i} = strcat(csv, ', ', ID(i).T);
    end 
    save(generateName(options, 'matfilein'), 'Spectra', 'CompleteSpectra', 'SpectraNames', '-append');

    fprintf('Reading MSI data according to ID file.\n');
    MSIs = cell(msiN,1);
    Masks = cell(msiN,1);
    MaskIs = cell(msiN,1);
    MSINames = cell(msiN,1);
    WhiteIs = cell(msiN,1);
    DarkIs = cell(msiN,1);

    groups = findgroups([ID.MsiID]);
    for g = 1:max(groups)
    %% Read MSI
        gIdxs = find(groups == g);
        gMembers = ID(gIdxs);
        coordinates = [gMembers.Originx; gMembers.Originy];

        currentOptions = options;
        names = cell(length(gMembers), 1);
        for j = 1:length(gMembers)
            [plotName, name] = generateName(options, 'read', ID(gIdxs(j)));
            currentOptions.saveOptions.plotName{j} = plotName;
            names{j} = name;
        end

        idd = ID(gIdxs(1));
        files = {data([data.MsiID] == idd.MsiID).File};    
        if contains(options.dataset, 'region')
            [segmentMSI, segmentWhite, segmentDark, segmentMask, segmentMaskI] = segmentMSIRegion(files, coordinates, currentOptions, 0.65, 35, 0.08);

        else % (square case)
            [segmentMSI, segmentWhite, segmentDark, segmentMask, segmentMaskI] = readMSI(files, coordinates, 5, 5, currentOptions); 
        end  

    %% Save MSI 
        for j = 1:length(gMembers)
            jj = gIdxs(j);
            MSIs{jj} = segmentMSI{j};
            Masks{jj} = segmentMask{j};
            MaskIs{jj} = segmentMaskI{j};
            MSINames{jj} = names{j};
            WhiteIs{jj} = segmentWhite{j};
            DarkIs{jj} = segmentDark{j};
        end
    end
    save(generateName(options, 'matfilein'), 'MSIs', 'Masks', 'MSINames', 'WhiteIs', 'DarkIs', '-append');
    save(generateName(options, 'matfilein-v7.3'), 'MaskIs', '-v7.3');        
    warning(w);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Precompute correlation smoothing Matrices 
    fprintf('Pre-computing smoothing matrices.\n');
    precomputedFile = fullfile(options.systemdir, 'precomputedParams.mat'); %pre-set parameters
    samples = unique({ID.Sample});

    %% Smoothing matrix from recorded spectra
    rAll = zeros(81, specN);
    rBenign = zeros(81, specN);
    rMalignant = zeros(81, specN);
    rFixed = zeros(81, specN);
    rCut = zeros(81, specN);
    rUnfixed = zeros(81, specN);
    rBenignCut = zeros(81, specN);
    rBenignFixed = zeros(81, specN);
    rBenignUnfixed = zeros(81, specN);
    rMalignantCut = zeros(81, specN);
    rMalignantFixed = zeros(81, specN);
    rMalignantUnfixed = zeros(81, specN);

    Cor_SampleBenignCut = zeros(length(samples), wavelengthN, wavelengthN);
    Cor_SampleBenignFixed = zeros(length(samples), wavelengthN, wavelengthN);
    Cor_SampleBenignUnfixed = zeros(length(samples), wavelengthN, wavelengthN);
    Cor_SampleMalignantCut = zeros(length(samples), wavelengthN, wavelengthN);
    Cor_SampleMalignantFixed = zeros(length(samples), wavelengthN, wavelengthN);
    Cor_SampleMalignantUnfixed = zeros(length(samples), wavelengthN, wavelengthN);
    Cor_Sample = zeros(length(samples), wavelengthN, wavelengthN);
    Cor_SampleMalignant = zeros(length(samples), wavelengthN, wavelengthN);
    Cor_SampleBenign = zeros(length(samples), wavelengthN, wavelengthN);

    total = 0;
    for i = 1:length(samples)
        %find the relative samples
        name = samples{i};
        sampleSpectra = uniqueSpectraNames(contains(uniqueSpectraNames, name));
        sampleSpectraIdx = find(contains(uniqueSpectraNames, name));

        sampleNum = numel(sampleSpectraIdx);
        rSample = zeros(81, sampleNum);
        rSampleBenign = zeros(81, sampleNum);
        rSampleMalignant = zeros(81, sampleNum);
        rSampleBenignCut = zeros(81, sampleNum);
        rSampleBenignFixed = zeros(81, sampleNum);
        rSampleBenignUnfixed = zeros(81, sampleNum);
        rSampleMalignantCut = zeros(81, sampleNum);
        rSampleMalignantFixed = zeros(81, sampleNum);
        rSampleMalignantUnfixed = zeros(81, sampleNum);

        for j = 1:length(sampleSpectra)
            total = total + 1;
            jj = sampleSpectraIdx(j);
            id = ID(uniqueSpectraIdxs(jj));
            ref = uniqueSpectra(jj, :)'; 

            rAll(:, total) =  ref;
            rSample(:, j) =  ref;

            if (id.IsBenign)
                %% Benign case
                rBenign(:, total) =  ref; 
                rSampleBenign(:, j) =  ref; 

                if (id.IsCut)
                    rCut(:, total) = ref; 
                    rBenignCut(:, total) = ref; 
                    rSampleBenignCut(:, j) = ref; 
                elseif (id.IsFixed)
                    rFixed(:, total) = ref; 
                    rBenignFixed(:, total) = ref; 
                    rSampleBenignFixed(:, j) = ref; 
                else
                    rUnfixed(:, total) = ref; 
                    rBenignUnfixed(:, total) = ref; 
                    rSampleBenignUnfixed(:, j) = ref; 
                end
            else
                %% Malignant case 
                rMalignant(:, total) = ref; 
                rSampleMalignant(:, j) = ref;

                if (id.IsCut)
                    rCut(:, total) = ref; 
                    rMalignantCut(:, total) = ref; 
                    rSampleMalignantCut(:, j) = ref; 
                elseif (id.IsFixed)
                    rFixed(:, total) = ref; 
                    rMalignantFixed(:, total) = ref;
                    rSampleMalignantFixed(:, j) = ref; 
                else
                    rUnfixed(:, total) = ref; 
                    rMalignantUnfixed(:, total) = ref; 
                    rSampleMalignantUnfixed(:, j) = ref; 
                end
            end
        end  

        Cor_SampleBenignCut(i,:,:) = smoothingMatrix(rSampleBenignCut); 
        Cor_SampleBenignFixed(i,:,:) = smoothingMatrix(rSampleBenignFixed); 
        Cor_SampleBenignUnfixed(i,:,:) = smoothingMatrix(rSampleBenignUnfixed);  

        Cor_SampleMalignantCut(i,:,:) = smoothingMatrix(rSampleMalignantCut);
        Cor_SampleMalignantFixed(i,:,:) = smoothingMatrix(rSampleMalignantFixed);
        Cor_SampleMalignantUnfixed(i,:,:) = smoothingMatrix(rSampleMalignantUnfixed);

        % correlation of all the measured spectra for this sample
        Cor_Sample(i,:,:) = smoothingMatrix(rSample); 

        % correlation of all cancerous spectra for this sample
        Cor_SampleMalignant(i,:,:) = smoothingMatrix(rSampleMalignant);

        % correlation of all normal spectra for this sample
        Cor_SampleBenign(i,:,:) = smoothingMatrix(rSampleBenign);

    end

    Cor_All = smoothingMatrix(rAll);
    Cor_Benign= smoothingMatrix(rBenign);
    Cor_Malignant = smoothingMatrix(rMalignant);
    Cor_Unfixed = smoothingMatrix(rUnfixed);
    Cor_Fixed = smoothingMatrix(rFixed);
    Cor_Cut = smoothingMatrix(rCut);
    Cor_BenignCut = smoothingMatrix(rBenignCut);
    Cor_BenignUnfixed = smoothingMatrix(rBenignUnfixed);
    Cor_BenignFixed = smoothingMatrix(rBenignFixed);
    Cor_MalignantCut = smoothingMatrix(rMalignantCut);
    Cor_MalignantUnfixed = smoothingMatrix(rMalignantUnfixed);
    Cor_MalignantFixed = smoothingMatrix(rMalignantFixed);

    save(precomputedFile, 'Cor_All', 'Cor_Benign', 'Cor_Malignant', 'Cor_Unfixed', 'Cor_Fixed',...
        'Cor_Cut', 'Cor_BenignCut', 'Cor_BenignUnfixed', 'Cor_BenignFixed', 'Cor_MalignantCut', ...
        'Cor_MalignantUnfixed', 'Cor_MalignantFixed', 'Cor_SampleBenignCut', 'Cor_SampleBenignFixed', ...
        'Cor_SampleBenignUnfixed', 'Cor_SampleMalignantCut', 'Cor_SampleMalignantFixed', ...
        'Cor_SampleMalignantUnfixed', 'Cor_Sample', 'Cor_SampleMalignant', 'Cor_SampleBenign', '-append');

    %% Count dataset
    datasetBreakdown(ID);

    %% Smoothing matrix based on markovian process 
    if isfield(options, 'rho')
        rho = options.rho;
    else 
        rho = 0.98;
    end

    M_row = ones(1,wavelengthN);
    for i=2:wavelengthN
        M_row(i) = M_row(i-1) * rho;
    end

    M = zeros(wavelengthN); %first order Markov process covariance matrix
    select = eye(1, wavelengthN);
    D = zeros(wavelengthN,wavelengthN); % Q x lambdaN
    for i=1:wavelengthN
        M(i,:) = [ fliplr(M_row(1: i)) , M_row(2:wavelengthN-i+1)];
        D(i,:) = circshift( select', (i-1)*wavelengthN); %diagonal with 1 in the position of components
    end
    Cor_Markovian = M;
    save(precomputedFile, 'Cor_Markovian', '-append');

    %% Smoothing matrix based on macbeth chart spectra
    load(precomputedFile, 'avgMeasuredMacbeth');
    macbeths = avgMeasuredMacbeth(wavelengthIdxs,:);

    z=xcorr(macbeths, macbeths, 'coeff');
    r = z(wavelengthN:end);
    Cor_Macbeth = toeplitz(r);
    save(precomputedFile, 'Cor_Macbeth', '-append');

    %% IlluminationXSensitivity
    Hgreen = illumXsensitivity(illumination, sensitivity, 'green');
    Hadjusted = illumXsensitivity(illumination, sensitivity, 'adjusted');
    Hrms = illumXsensitivity(illumination, sensitivity, 'rms');
    Hextended = illumXsensitivity(illumination, sensitivity, 'extended');
    save(precomputedFile, 'Hgreen', 'Hadjusted', 'Hrms', 'Hextended', '-append');

    %% Create coefficient table
    disp('Reading coefficients from excel file...')
    pixelValueSelectionMethods = {'green', 'rms', 'adjusted'};
    Coefficients = ones(length(ID), 3, 7);
    xls = dir(fullfile(options.systemdir, 'coeff*'));
    src = fullfile(options.systemdir, xls(1).name);
    for k = 1:msiN
        % Retrieve MSI data
        g = MSIs{k};
        % Retrieve spectrum data
        refmeasured = CompleteSpectra(uniqueSpectraIdxsInID(k),:); 

        for m = 1:length(pixelValueSelectionMethods)
            gg = raw2msi(g, pixelValueSelectionMethods{m});
            [r, c] = find(MaskIs{k});
            x = ID(k).Originx - min(c) + 1;
            y = ID(k).Originy - min(r) + 1;
            refmsi = im2uint16(squeeze(gg(:,y,x)))';
            xlswrite2(src, refmeasured, 4, 'B3');
            xlswrite2(src, refmsi, 4, 'D407:J407');
            ctemp = xlsread(src, 4, 'D414:J414');
            Coefficients(k, m, :) = ctemp;
        end
    end
    save(fullfile(options.systemdir, 'precomputedParams.mat'), 'Coefficients', '-append');
    load(fullfile(options.systemdir, 'precomputedParams.mat'), 'Coefficients');
    disp('Searching for reference coefficient...')
    ID = selectCoeffIndex(options); 
else   
    %% Read newly created files 

    load(generateName(options, 'matfilein'));
    name = options.dataset;
end
    

fprintf('Finished initialization.\n\n')
fprintf('Let''s get down to business :muscle:\n\n');
% LOAD DATA & INITIALIZE ends   
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Support function for computing the smoothing matrix from counts
function M = smoothingMatrix(s)
%smoothingMatrix = @(s) 1 / (size(s, 2) - 1) * s * s';  %@(s) 1 / (size(s, 2) - 1) * (s -  mean(s,2)) * (s -  mean(s,2))';
    sTrue = s(:,~all(s == 0, 1));
    M = 1 / size(sTrue, 2) * (sTrue * sTrue');
end

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