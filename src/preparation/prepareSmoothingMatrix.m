
%% make correlation smoothing matrices

in = matfile( fullfile(options.systemdir, 'in.mat')); 
param = matfile(fullfile(options.systemdir, 'precomputedParams.mat'), 'Writable', true);

samples = unique([ID.Sample]);
[spectra, spectraIdx, ~] = unique(strcat({ID.Csvid}, {ID.T}));
dataNum = numel(spectraIdx);

wavelengthN = size(sensitivity, 1);
wavelength = linspace(380, 780, wavelengthN);

%% Smoothing matrix from recorded spectra
rAll = zeros(81, dataNum);
rBenign = zeros(81, dataNum);
rMalignant = zeros(81, dataNum);
rFixed = zeros(81, dataNum);
rCut = zeros(81, dataNum);
rUnfixed = zeros(81, dataNum);
rBenignCut = zeros(81, dataNum);
rBenignFixed = zeros(81, dataNum);
rBenignUnfixed = zeros(81, dataNum);
rMalignantCut = zeros(81, dataNum);
rMalignantFixed = zeros(81, dataNum);
rMalignantUnfixed = zeros(81, dataNum);

for i = 1:length(samples)
    %find the relative samples
    name = samples{i};
    sampleSpectra = spectra(contains(spectra, name));
    sampleSpectraIdx = spectraIdx(contains(spectra, name));
    
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
      
        id = ID(sampleSpectraIdx(j));
        
        ref = interp1(380:780, measuredSpectrumStruct(j).Spectrum, wavelength, 'nearest')';
        rAll = [rAll, ref];
        rSample = [rSample, ref];
       
        if (id.IsNormal)
            %% Benign case
            rBenign = [rBenign, ref];
            rSampleBenign = [rSampleBenign, ref];
            
            if (id.IsCut)
                rCut = [rCut, ref];
                rBenignCut = [rBenignCut, ref];
                rSampleBenignCut = [rSampleBenignCut, ref];
            elseif (id.IsFixed)
                rFixed = [rFixed, ref];
                rBenignFixed = [rBenignFixed, ref];
                rSampleBenignFixed = [rSampleBenignFixed, ref];
            else
                rUnfixed = [rUnfixed, ref];
                rBenignUnfixed = [rBenignUnfixed, ref];
                rSampleBenignUnfixed = [rSampleBenignUnfixed, ref];
            end
        else
            %% Malignant case 
            rMalignant = [rMalignant, ref];
            rSampleMalignant = [rSampleMalignant, ref];
            
            if (id.IsCut)
                rCut = [rCut, ref];
                rMalignantCut = [rMalignantCut, ref];
                rSampleMalignantCut = [rSampleMalignantCut, ref];
            elseif (id.IsFixed)
                rFixed = [rFixed, ref];
                rMalignantFixed = [rMalignantFixed, ref];
                rSampleMalignantFixed = [rSampleMalignantFixed, ref];
            else
                rUnfixed = [rUnfixed, ref];
                rMalignantUnfixed = [rMalignantUnfixed, ref];
                rSampleMalignantUnfixed = [rSampleMalignantUnfixed, ref];
            end
        end
    end  
    
    param.(matlab.lang.makeValidName(['Cor_SampleBenignCut_', name])) = smoothingMatrix(rSampleBenignCut); 
    param.(matlab.lang.makeValidName(['Cor_SampleBenignFixed_', name])) = smoothingMatrix(rSampleBenignFixed); 
    param.(matlab.lang.makeValidName(['Cor_SampleBenignUnfixed_', name])) = smoothingMatrix(rSampleBenignUnfixed);  
    
    param.(matlab.lang.makeValidName(['Cor_SampleMalignantCut_', name])) = smoothingMatrix(rSampleMalignantCut);
    param.(matlab.lang.makeValidName(['Cor_SampleMalignantFixed_', name])) = smoothingMatrix(rSampleMalignantFixed);
    param.(matlab.lang.makeValidName(['Cor_SampleMalignantUnfixed_', name])) = smoothingMatrix(rSampleMalignantUnfixed);
    
    % correlation of all the measured spectra for this sample
    param.(matlab.lang.makeValidName(['Cor_', name])) = smoothingMatrix(rSample); 
    
    % correlation of all cancerous spectra for this sample
    param.(matlab.lang.makeValidName(['Cor_SampleMalignant_', name])) = smoothingMatrix(rSampleMalignant);
    
    % correlation of all normal spectra for this sample
    param.(matlab.lang.makeValidName(['Cor_SampleBenign_', name])) = smoothingMatrix(rSampleBenign);

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

%% Support function for computing the smoothing matrix from counts
function K = smoothingMatrix(s)
    means = mean(s,2);
    K = 1 / (size(s, 2) - 1) * (s - means) * (s - means)';
% %     alternatively 
%     z = xcorr(s, s, 'coeff');
%     r = z(wavelengthN:end);
%     K = toeplitz(r);
end
    
