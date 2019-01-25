%% LOAD DATA & INITIALIZE
if ~options.skipLoading
    
    disp('Initialization.')
    load(fullfile(options.systemdir, 'system.mat')); % camera system parameters
    load(fullfile(options.systemdir, 'data.mat')); % image data
    load(fullfile(options.systemdir, 'ID.mat')); % image data id and info struct

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    msiN = length(ID);
    if (options.tryReadData)
        
        w = warning('off', 'all');
        
        fprintf('Reading spectral data according to ID file.\n');
        m = matfile(generateName(options, 'matfilein'), 'Writable', true);

        %% Read raw spectra 
        [uniqueSpectraNames, uniqueSpectraIdxs, uniqueSpectraIdxsInID] = unique(strcat({ID.Csvid}, {ID.T}));
        specN = length(uniqueSpectraIdxs);
        uniqueSpectra = zeros(401, specN);

        for i = 1:specN
            % read raw measured spectrum
            rawSpectrum = readSpectrum(ID(i).Csvid, ID(i).T);
            % read raw white measured spectrum of the reference surface
            referenceSpectrum = readSpectrum(char(strcat(data(ID(i).Representative).Sample, '\', 'white.csv')));

            if abs(rawSpectrum-referenceSpectrum) < 0.000001
                error('Measurement is same as white.')
            end

            uniqueSpectra(:,i) = rawSpectrum ./ referenceSpectrum;

        end
         
        for i = 1:msiN
            
        %% Save Spectra Struct
            if ~isempty(uniqueSpectra(:, uniqueSpectraIdxsInID))
                [~, csv] =  generateName(options, 'csv', ID(i));
                m.MeasuredSpectrumStruct(i,1) = struct( 'Index', i, 'Name', csv, 'Spectrum', uniqueSpectra(:, uniqueSpectraIdxsInID(i)), 'T', ID(i).T);
            end
        end 
        
        fprintf('Reading MSI data according to ID file.\n');
        groups = findgroups([ID.UniqueCount]);
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
                directory = fileparts(plotName);
                if ~exist(directory, 'dir')
                    mkdir(directory);
                    addpath(directory);
                end
            end
            
            idd = ID(gIdxs(1));
            files = {data(idd.Data).File};    
            if contains(options.dataset, 'region')
                segment = segmentMSIRegion(files, coordinates, currentOptions);

            else % (square case)
                segment = readMSI(files, coordinates, 5, 5, currentOptions); 
            end  

        %% Save MSI   
            for j = 1:length(gMembers)
                jj = gIdxs(j);
                m.MSIStruct(jj,1) = struct('Name', names{j}, 'Index', jj, 'MSI', segment(j).MSI, 'Mask', segment(j).patchMask, 'MaskI', segment(j).maskI);
                m.WhiteMSIStruct(jj,1) = struct('Name', names{j}, 'Index', jj, 'MSI', segment(j).whiteReference);
                m.DarkMSIStruct(jj,1) = struct('Name', names{j}, 'Index', jj, 'MSI', segment(j).darkReference);
            end
        end
        
        warning(w);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Precompute correlation smoothing Matrices 
        fprintf('Pre-computing smoothing matrices.\n');
        param = matfile(fullfile(options.systemdir, 'precomputedParams.mat'), 'Writable', true);

        samples = unique([ID.Sample]);
        wavelengthN = size(sensitivity, 1);
        wavelength = linspace(380, 780, wavelengthN);

        %% Support function for computing the smoothing matrix from counts
        smoothingMatrix = @(s) 1 / (size(s, 2) - 1) * (s -  mean(s,2)) * (s -  mean(s,2))';

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
                ref = interp1(380:780, uniqueSpectra(:, jj), wavelength, 'nearest')';

                rAll(:, total) =  ref; %[rAll, ref];
                rSample(:, j) =  ref; %[rSample, ref];

                if (id.IsNormal)
                    %% Benign case
                    rBenign(:, total) =  ref; %[rBenign, ref];
                    rSampleBenign(:, j) =  ref; %[rSampleBenign, ref];

                    if (id.IsCut)
                        rCut(:, total) = ref; %[rCut, ref];
                        rBenignCut(:, total) = ref; %[rBenignCut, ref];
                        rSampleBenignCut(:, j) = ref; %[rSampleBenignCut, ref];
                    elseif (id.IsFixed)
                        rFixed(:, total) = ref; %[rFixed, ref];
                        rBenignFixed(:, total) = ref; %[rBenignFixed, ref];
                        rSampleBenignFixed(:, j) = ref; %[rSampleBenignFixed, ref];
                    else
                        rUnfixed(:, total) = ref; %[rUnfixed, ref];
                        rBenignUnfixed(:, total) = ref; %[rBenignUnfixed, ref];
                        rSampleBenignUnfixed(:, j) = ref; %[rSampleBenignUnfixed, ref];
                    end
                else
                    %% Malignant case 
                    rMalignant(:, total) = ref; %[rMalignant, ref];
                    rSampleMalignant(:, j) = ref; %[rSampleMalignant, ref];

                    if (id.IsCut)
                        rCut(:, total) = ref; %[rCut, ref];
                        rMalignantCut(:, total) = ref; %[rMalignantCut, ref];
                        rSampleMalignantCut(:, j) = ref; %[rSampleMalignantCut, ref];
                    elseif (id.IsFixed)
                        rFixed(:, total) = ref; %[rFixed, ref];
                        rMalignantFixed(:, total) = ref; %[rMalignantFixed, ref];
                        rSampleMalignantFixed(:, j) = ref; %[rSampleMalignantFixed, ref];
                    else
                        rUnfixed(:, total) = ref; %[rUnfixed, ref];
                        rMalignantUnfixed(:, total) = ref; %[rMalignantUnfixed, ref];
                        rSampleMalignantUnfixed(:, j) = ref; %[rSampleMalignantUnfixed, ref];
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

        param.(matlab.lang.makeValidName('Cor_All')) = smoothingMatrix(rAll);
        param.(matlab.lang.makeValidName('Cor_Benign')) = smoothingMatrix(rBenign);
        param.(matlab.lang.makeValidName('Cor_Malignant')) = smoothingMatrix(rMalignant);
        param.(matlab.lang.makeValidName('Cor_Unfixed')) = smoothingMatrix(rUnfixed);
        param.(matlab.lang.makeValidName('Cor_Fixed')) = smoothingMatrix(rFixed);
        param.(matlab.lang.makeValidName('Cor_Cut')) = smoothingMatrix(rCut);
        param.(matlab.lang.makeValidName('Cor_BenignCut')) = smoothingMatrix(rBenignCut);
        param.(matlab.lang.makeValidName('Cor_BenignUnfixed')) = smoothingMatrix(rBenignUnfixed);
        param.(matlab.lang.makeValidName('Cor_BenignFixed')) = smoothingMatrix(rBenignFixed);
        param.(matlab.lang.makeValidName('Cor_MalignantCut')) = smoothingMatrix(rMalignantCut);
        param.(matlab.lang.makeValidName('Cor_MalignantUnfixed')) = smoothingMatrix(rMalignantUnfixed);
        param.(matlab.lang.makeValidName('Cor_MalignantFixed')) = smoothingMatrix(rMalignantFixed);

        param.(matlab.lang.makeValidName('Cor_SampleBenignCut')) = Cor_SampleBenignCut;
        param.(matlab.lang.makeValidName('Cor_SampleBenignFixed')) = Cor_SampleBenignFixed;
        param.(matlab.lang.makeValidName('Cor_SampleBenignUnfixed')) = Cor_SampleBenignUnfixed;
        param.(matlab.lang.makeValidName('Cor_SampleMalignantCut')) = Cor_SampleMalignantCut;
        param.(matlab.lang.makeValidName('Cor_SampleMalignantFixed')) = Cor_SampleMalignantFixed;
        param.(matlab.lang.makeValidName('Cor_SampleMalignantUnfixed')) = Cor_SampleMalignantUnfixed;
        param.(matlab.lang.makeValidName('Cor_Sample')) = Cor_Sample;
        param.(matlab.lang.makeValidName('Cor_SampleMalignant')) = Cor_SampleMalignant;
        param.(matlab.lang.makeValidName('Cor_SampleBenign')) = Cor_SampleBenign;

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
        param.(matlab.lang.makeValidName('Cor_Markovian')) = M;

        %% Smoothing matrix based on macbeth chart spectra
        macbeths = interp1(380:780, param.avgMeasuredMacbeth, wavelength, 'nearest')';

        z=xcorr(macbeths, macbeths, 'coeff');
        r = z(wavelengthN:end);
        param.(matlab.lang.makeValidName('Cor_Macbeth')) = toeplitz(r);

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Read newly created files 
    in = matfile(generateName(options, 'matfilein')');
    MSIStruct = in.MSIStruct;
    whiteStruct = in.WhiteMSIStruct;
    darkStruct = in.DarkMSIStruct;
    measuredSpectrumStruct = in.MeasuredSpectrumStruct;
    name = options.dataset;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Create coefficient table
    if (options.tryReadData)
        pixelValueSelectionMethods = {'green', 'rms', 'adjusted'};
        coeff = ones(length(ID), 3, 7);
        src = '../../input/others/coeff.xlsx';
        for k = 1:msiN
            % Retrieve MSI data
            g = MSIStruct(k).MSI;

            % Retrieve spectrum data
            refmeasured = measuredSpectrumStruct(k).Spectrum;

            for m = 1:length(pixelValueSelectionMethods)
                gg = valueSelect(g, pixelValueSelectionMethods{m});
                if contains(options.systemdir, 'region')
                    mask = reshape( MSIStruct(k).Mask, 1, size(gg,2)*size(gg,3));
                    gg = reshape(gg, size(gg,1), size(gg,2)*size(gg,3));
                    gg = gg(:,mask);
                else
                    gg = gg(:,2:4,2:4); %optional step
                    gg = reshape(gg, size(gg,1), 3*3);
                end
                refmsi = im2uint16(mean(gg, 2))';
                xlswrite2(src, refmeasured, 4, 'B3');
                xlswrite2(src, refmsi, 4, 'D407:J407');
                ctemp = xlsread(src, 4, 'D414:J414');
                coeff(k, m, :) = ctemp;
            end
        end
        save(fullfile(options.systemdir, 'coeff.mat'), 'coeff', '-v7.3');
        ID = selectCoeffIndex(options);  
    end
    
end

fprintf('Finished initialization.\n\n')
fprintf('Let''s get down to business :muscle:\n\n');
% LOAD DATA & INITIALIZE ends   
    
