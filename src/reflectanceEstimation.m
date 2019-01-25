function [estimatedReflectance, rmse, nmse, minIdx] = reflectanceEstimation(g, spectrum, id, options)

%% wienerEstimation Performs wiener estimation for reflectance from images and returns an sRGB reconstucted image
%Input arguments
%g = the input MSI struct (name, index, 4D-MSI, RegionMask), where 4D-MSI(filter band, x , y , image rgb values)
%spectrum = the measured originalReflectance spectrum of the sample
%id = the point id information
%
%Optional input arguments
%'pixelValueSelectionMethod' = {'green', 'rms', 'adjusted', 'extended', 'unchanged'} for selecting a single value from the 3 RGB values
%'smoothingMatrixMethod' = {'markovian 0.99', 'adaptive', 'KCor all data', 'KCor macbeth', 'KCor all specimen', ...
%                'KCor all same malignancy', 'KCor all same fixation' 'KCor same malignancy', 'KCor same malignancy, fixation'}; for selecting the smoothing matrix
%'rho' =  parameter for markovian smoothing matrix , default: 0.99
%'variance' =  parameter for noise covariance matrix , default: 1
%'noiseType' = to include noise component,{'none', 'sameForChannel', 'difForChannel', 'spatial', 'givenSNR', 'white gaussian', 'fromOlympus'}
%'reference' = the reference white image values
%'windowDim' = the dimension of square pixel neighborhood (odd number), default: 3
%'SVDTol' = to use tolerance for SVD when computing an inverse matrix, default: false 
%
%Output arguments
%estimatedReflectance = the estimated reflectance for every pixel in the
%sample area
%meanEstimatedReflectance = the average normalized reflectance estimation of the sample area
%rmse = root mean square error
%nmse = normalized mean square error

%% Argument parsing

c = matfile(fullfile(options.systemdir, 'system.mat'));
sensitivity = c.sensitivity; %sensitivity = spectral sensitivity of RGB camera (3 columns)
illumination = c.illumination; %illumination = spectral illumination for all channels
clear('c');

if isstruct(g)
    MSI = g.MSI; 
    if ~isfield(g, 'Mask')
        mask = ones(size(MSI, 2), size(MSI, 3));
    else
        mask = g.Mask;
    end
else 
    MSI = g;
    mask = ones(size(g, 2), size(g, 3));
end

m = matfile(fullfile(options.systemdir, 'precomputedParams.mat')); %pre-set parameters

wavelengthN = size(sensitivity, 1);
wavelength = linspace(380, 780, wavelengthN);

if isfield(options, 'smoothingMatrixMethod')
    smoothingMatrixMethod = options.smoothingMatrixMethod;
else 
    smoothingMatrixMethod = 'markovian';
end

if strcmp(smoothingMatrixMethod, 'markovian')
    defaultRho = 0.99;
    if isfield(options, 'rho')
        rho = options.rho;
    else 
        rho = defaultRho;
    end
end

if strcmp(smoothingMatrixMethod, 'adaptive')
    defaultAlpha = 0.4;
    if isfield(options, 'alpha')
        alpha = options.alpha;
    else
        alpha = defaultAlpha;
    end
end 

if isfield(options, 'pixelValueSelectionMethod')
    pixelValueSelectionMethod = options.pixelValueSelectionMethod;
else 
    pixelValueSelectionMethod = 'extended';
end

if isfield(options, 'noiseType')
    noiseType = options.noiseType;
else 
    noiseType = 'independent';
end

if strcmp(pixelValueSelectionMethod, 'rgb')
    illumination = illumination(:, 8); % white light
else
    illumination = illumination(:, 1:7); % narrow bandwidth lights
end


if strcmp(options.pixelValueSelectionMethod, 'rgb')
    coeff = ones(3, 1) / 10^3;
else
    c = matfile(fullfile(options.systemdir, 'coeff.mat'));
    pixelValueSelectionMethods = {'green', 'rms', 'adjusted', 'extended', 'rgb'};
    coeff = squeeze(c.coeff(id.CoeffIndex, min(find(strcmp(pixelValueSelectionMethods, options.pixelValueSelectionMethod)), 3), :))';
    clear('c');
end

if strcmp(pixelValueSelectionMethod, 'extended')
    if ~isempty(coeff)
        coeff = [coeff(1:2), coeff(2), coeff(3:5), coeff(6), coeff(6:7)];
    end
end

if isfield(options, 'SVDTol')
    hasSVDTol = true;
    tol = options.SVDTol;
else
    hasSVDTol = false;
    tol = 1;
end

% sample = generateName([], 'sample', [], id);

isBenign = id.IsNormal;
if (id.IsCut)
    fixing = 'Cut';
elseif (id.IsFixed)
    fixing = 'Fixed';
else
    fixing = 'Unfixed';
end

% Argument parsing ends

%% Generate smoothing matrix for estimation

G = valueSelect(MSI, pixelValueSelectionMethod); % convert from 4D MSI+rgb to 3D MSI+grey
[msibands, height, width] = size(G);
G = G ./ options.luminanceCorrection;

coeff = coeff * 5; % to adapt coefficients from bandwidth 1nm to bandwidth 5nm (Source is 1x141, esimation is 1x81)

% computed beforehand with prepareSmoothingMatrix
switch smoothingMatrixMethod
    case {'Markovian 0.99', 'markovian'}
        if rho ~= defaultRho
            error('Invalid correlation parameter for the markovian process');
        end
        M = m.Cor_Markovian;
        M = mean(m.signalCov) * M;
        
    case 'Cor_All' % average xcorr of all measured data
        M = m.Cor_All;
        
    case 'Cor_Macbeth' % average xcorr of all macbeth data
        M = m.Cor_Macbeth;
        
    case 'Cor_Sample' % average xcorr of all measured data for each sample
        M = squeeze(m.Cor_Sample(id.SampleId,:,:));
        
    case 'Cor_Malignancy'
        if isBenign
            M = m.Cor_Benign;  %m.Mn; 
            
        else
            M = m.Cor_Malignant; %m.Mc;
        end
        
    case 'Cor_Fixing'
        if strcmp(fixing, 'Cut')
            M = m.Cor_Cut;
        elseif strcmp(fixing, 'Fixed')
            M = m.Cor_Fixed;
        else
            M = m.Cor_Unfixed;
        end
        
    case 'Cor_SampleMalignancy' % average xcorr of all measured cancerous data for cancerous sample
        if isBenign
            M = squeeze(m.Cor_SampleBenign(id.SampleId,:,:));
        else
            M = squeeze(m.Cor_SampleMalignant(id.SampleId,:,:));
        end
        
    case 'Cor_SampleMalignancyFixing' % average xcorr of all fixed cancer measured data for fixed cancer sample
        if isBenign && strcmp(fixing, 'Cut')
            M = squeeze(m.Cor_SampleBenignCut(id.SampleId,:,:));
        elseif isBenign && strcmp(fixing, 'Fixed')
            M = squeeze(m.Cor_SampleBenignFixed(id.SampleId,:,:));
        elseif isBenign && strcmp(fixing, 'Unfixed')
            M = squeeze(m.Cor_SampleBenignUnfixed(id.SampleId,:,:));
        elseif ~isBenign && strcmp(fixing, 'Cut')
            M = squeeze(m.Cor_SampleMalignantCut(id.SampleId,:,:));
        elseif ~isBenign && strcmp(fixing, 'Fixed')
            M = squeeze(m.Cor_SampleMalignantFixed(id.SampleId,:,:));
        elseif ~isBenign && strcmp(fixing, 'Unfixed')
            M = squeeze(m.Cor_SampleMalignantUnfixed(id.SampleId,:,:));
        else 
            error('Unsupported type')
        end
        
    case 'Cor_MalignancyFixing'
        if isBenign && strcmp(fixing, 'Cut')
            M = m.Cor_BenignCut;
        elseif isBenign && strcmp(fixing, 'Fixed')
            M = m.Cor_BenignFixed;
        elseif isBenign && strcmp(fixing, 'Unfixed')
            M = m.Cor_BenignUnfixed;
        elseif ~isBenign && strcmp(fixing, 'Cut')
            M = m.Cor_MalignantCut;
        elseif ~isBenign && strcmp(fixing, 'Fixed')
            M = m.Cor_MalignantFixed;
        elseif ~isBenign && strcmp(fixing, 'Unfixed')
            M = m.Cor_MalignantUnfixed;
        else 
            error('Unsupported type')
        end
        
    case 'adaptive'
        %% Based on "Reflectance reconstruction for multispectral imaging by adaptive Wiener estimation"[Shen2007]
        options.smoothingMatrixMethod = 'Cor_Sample';
        rhat = reflectanceEstimation(MSI, spectrum, id, options);
        M = adaptiveSmoothingMatrix(rhat, options.systemdir, alpha);
        
    case 'spatiospectral'
%         %unfinished      
        
    case 'obi2001'
%         %unfinished      

    otherwise
        error('Unexpected smoothing matrix method. Abort execution.')
end

% Generate smoothing matrix for estimation ends

%% Illumination X sensitivity

H = illumXsensitivity(illumination, sensitivity, pixelValueSelectionMethod); % illumination * sensitivity
H = diag(coeff) * H;

MH = M * H';
HMH = H * M * H';

% Illumination X sensitivity ends

%% Covariance matrix of the additive noise
if contains(noiseType, 'sameForChannel')
    attr = strsplit(noiseType, {' '});
    if (numel(attr) > 1)
        sigma = str2double(attr{2});
    elseif isfield(options, 'sigma')
        sigma = options.sigma;
    else
        defaultSigma = 0.0001;
        sigma = defaultSigma;
    end
    variance = sigma *ones(msibands, 1);
    
elseif contains(noiseType, 'diffForChannel')
    attr = strsplit(noiseType, {' ', '[', ',', ']'});
    if (numel(attr) > 1)
        variance = (cellfun(@str2num,  attr(1, 2:(end-1)))).^2;
    elseif isfield(options, 'sigma')
        variance = options.sigma.^2;
    else
        variance = [0.0031, 0.0033, 0.0030, 0.0031, 0.0032, 0.0029, 0.0024];
    end
    
elseif contains(noiseType, 'givenSNR')
    attr = strsplit(noiseType, {' ', 'dB'});
    if (numel(attr) > 1)
        snr = str2double(attr{2});
    elseif isfield(options, 'snr')
        snr = options.snr;
    else
        defaultSNR = 17.5; %dB
        snr = defaultSNR; %dB
    end
    variance = (trace(HMH) / (msibands * 10^(snr / 10))) * ones(msibands, 1);
    
elseif contains(noiseType, 'white gaussian')
    attr = strsplit(options.noiseType, {' ', '^{', '}'});
    if (numel(attr) > 1)
        options.noiseOrder = str2double(attr{3})^str2double(attr{4});
    end
    variance = (randn(1, msibands) .* options.noiseOrder).^2;  %different at every channel

elseif strcmp(noiseType, 'fromOlympus')
    variance = (m.noiseparam').^2;
    
elseif strcmp(noiseType, 'none')
    variance = zeros(1, msibands);
    
elseif strcmp(noiseType, 'spatial') || contains(noiseType, 'spatial')
    if isfield(options, 'windowDim')
        windowDim = options.windowDim;
    else 
        windowDim = 5;
    end
    windowKernel = ones(windowDim);
    windowElements = windowDim^2;
    
    if size(windowKernel, 1) > height || size(windowKernel,2) > width
        windowKernel = ones(min(height, width));
        windowElements = min(height, width);
    end
    
    attr = strsplit(options.noiseType, {' '});
    
    if isfield(options, 'sigma1')
        sigma1 = options.sigma1;
    elseif (numel(attr) > 1)
        sigma1 = str2double(attr{2});
    end
    
    if isfield(options, 'sigma2')
        sigma2 = options.sigma2;
    elseif (numel(attr) > 2)
        sigma2 = str2double(attr{3});
    end 
    
    if exist('sigma1', 'var') && exist('sigma2', 'var')
        variance = ones(msibands,1) * (sqrt(0.5) * sigma1 + sigma2)^2; %else use noise variance from Olympus 
    elseif isfield(options, 'sigma') 
        variance = options.sigma.^2;
    else 
        variance = (m.noiseparam').^2;
    end
        
else 
    error('Not implemented');
end

if (length(variance) < msibands)
    variance = [variance(1:2), variance(2), variance(3:5), variance(6), variance(6:7)];
end
Kn = diag(variance); %different at every channel
    
% Covariance matrix of the additive noise ends

Gres = reshape(G, msibands, height*width); % msi = im2double(diag(coeff) * double(im2uint16(G(:,row,col))) );
mask = reshape(mask, 1, height*width);
    
if strcmp(noiseType, 'spatial')
%% Perform Spatially Adaptive Wiener estimation for all pixels in an image area
%  Based on "A Spatially AdaptiveWiener Filter for Reflectance Estimation"[Urban2008]

    meanG = zeros(msibands, width * height);
    centeredG = zeros(msibands, width * height);
    Kcov = zeros(msibands, width * height);
    for i = 1:msibands
        Gi = squeeze(G(i,:,:));
        means = conv2(Gi, windowKernel ./ windowElements,'same');        
        meanG(i, :) = reshape(means, 1, width * height);
        centeredG(i, :) =  reshape(Gi - means, 1, width * height);
        sdGi = stdfilt(Gi);
        varGi = sdGi.^2;
        Kcov(i, :) =  reshape(conv2(varGi, windowKernel ./ (windowElements - 1), 'same'), 1, width * height); 
    end
    
    activeRegion = find(mask);
    estimatedReflectance = zeros(wavelengthN, length(activeRegion));
    for i = 1:activeRegion
        meanGij = squeeze(meanG(:,i));
        centeredGij = squeeze(centeredG(:,i));
        Kij = diag(squeeze(Kcov(:,i)));
        W = multiplyToInverse(Kij, Kij + Kn); % same as kInWindow * inv(kInWindow + Kn)
        What = multiplyToInverse(MH, HMH + Kij - W * Kij, hasSVDTol, tol);
        estimatedReflectance(:, i) = What * W * centeredGij + What * meanGij;  
    end
    
else    
    activeRegion = Gres(:,mask);
    div = multiplyToInverse(MH , HMH + Kn, hasSVDTol, tol); % M 401x401, H' 401x7, inv() 7x7
    estimatedReflectance = div * activeRegion; % 401 x 100 
end  


%% Perform Wiener estimation for all pixels in an image area

if (height > 200 ||  width > 200) %in this case the mask is ones(height, width)
    estimatedReflectance = max(estimatedReflectance, 0);
    estimatedReflectance = min(estimatedReflectance, 1);
    
    estimatedReflectance = reshape(estimatedReflectance, numel(wavelength), height, width);
else   

    if any(estimatedReflectance(:) < 0 ) || any(estimatedReflectance(:) > 1 )
        warning('Estimated reflectance spectrum is out of bounds!');
    end
    
    idx = any(estimatedReflectance < 0) | any(estimatedReflectance > 1);
    if sum(idx) == size(estimatedReflectance, 2)
        estimatedReflectance = max(estimatedReflectance, 0);
        estimatedReflectance = min(estimatedReflectance, 1);
        warning('Estimated reflectance spectrum is out of bounds for all pixels in the region!')
    else
        estimatedReflectance = estimatedReflectance(:,~idx);
    end
    
    if ~isempty(spectrum)    
        [rmse, rmseIdx] = min(Rmse(spectrum, estimatedReflectance));
        [nmse, nmseIdx] = min(Nmse(spectrum, estimatedReflectance));
        %either min error index produces similar results. 
        estimatedReflectance = estimatedReflectance(:,rmseIdx);
        mini = activeRegion(rmseIdx) / width;
        minj = mod(activeRegion(rmseIdx) , width);
        minIdx = [mini, minj];
    end
end

% Perform the estimation for all pixels in an image area ends

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rmse = Rmse(r_measured, r_estimated)
        % Root Mean Square Error
        [N, M] = size(r_estimated);
        rmse = zeros(1, M);
        for i = 1:M
            rmse(i) = sqrt((r_measured - r_estimated(:,i))'*(r_measured - r_estimated(:,i))/ N);
        end
end

function nmse = Nmse(r_measured, r_estimated)
        % Normalized Mean Square Error
        [~, M] = size(r_estimated);
        nmse = zeros(1, M);
        for i = 1:M
            nmse(i) = (r_measured - r_estimated(:,i))' * (r_measured - r_estimated(:,i)) / (r_measured' * r_measured);
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subfunction for handling matrix inversion and multiplication
function C = multiplyToInverse(A, B, hasSVDTol, tol)

    if (nargin < 3)
        hasSVDTol = false;
    end

    if hasSVDTol 
        C = A * pinv(B, tol); %SVD using "tol" as singular value threshold
    else
        C = A / B; %same as A*inv(B)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function di = reflectanceDistance(ri, rhat, a)
    di = a * mean(abs(ri / norm(ri) - rhat / norm(rhat)))...
        + (1-a) * max(abs(ri / norm(ri) - rhat / norm(rhat)));
end

function k = replicationTimes(d, dmax, gamma)
    if (nargin < 3)
        gamma = 1;
    end
    k = floor((dmax / d)^gamma + 0.5);
end

function adaptedM = adaptiveSmoothingMatrix(rhat, systemdir, a)

    ff = matfile(fullfile(systemdir, 'in.mat'));
    measuredSpectra = ff.MeasuredSpectrumStruct;
    [~, idxs] = unique(strcat({measuredSpectra.Name}, {measuredSpectra.T}));
    measuredSpectra = measuredSpectra(idxs);
    wavelength = linspace(380, 780, 81)';
    n = length(measuredSpectra);
    r = zeros(81,n);
    d = zeros(1,n);
    for i = 1:n
        ri =  interp1(380:780, measuredSpectra(i).Spectrum, wavelength, 'nearest');
        d(i) = reflectanceDistance(ri, rhat, a);
        r(:,i) = ri;
    end
    dmax = max(d);
    
    spectra = [];
    for i = 1:n 
        if dmax > 0 
            k = replicationTimes(d(i), dmax);
        else 
            k = 1;
        end
        spectra = [ spectra, repmat(r(:,i),1,k)];
    end
    
    means = mean(spectra,2);
    adaptedM = 1 / (size(spectra, 2) - 1) * (spectra - means) * (spectra - means)';
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

    if (size(E, 1) ~= size(S, 1))

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