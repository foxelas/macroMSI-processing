function [] = computeSmoothingMatrices(measuredSpectra, measuredSpectraNames, wavelengthN, ID, options)
    %% Precompute correlation smoothing Matrices 
    fprintf('Pre-computing smoothing matrices.\n');
    precomputedFile = fullfile(options.systemdir, 'precomputedParams.mat'); %pre-set parameters
    samples = unique({ID.Sample});

    specN = length(measuredSpectraNames);
    wavelengthStep = ceil(length(380:780) / wavelengthN);
    wavelengthIdxs = 1:wavelengthStep:length(380:780);
    
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
        sampleSpectra = measuredSpectraNames(contains(measuredSpectraNames, name));
        sampleSpectraIdx = find(contains(measuredSpectraNames, name));

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
            id = ID(j);
            ref = measuredSpectra(jj, :)'; 

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
    breakdownDataset(ID);

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
end

%% Support function for computing the smoothing matrix from counts
function M = smoothingMatrix(s)
%smoothingMatrix = @(s) 1 / (size(s, 2) - 1) * s * s';  %@(s) 1 / (size(s, 2) - 1) * (s -  mean(s,2)) * (s -  mean(s,2))';
    sTrue = s(:,~all(s == 0, 1));
    M = 1 / size(sTrue, 2) * (sTrue * sTrue');
end
