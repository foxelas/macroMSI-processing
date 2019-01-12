function [estimatedReflectance, rmse, nmse] = reflectanceEstimation(g, spectrum, id, options)

%% wienerEstimation Performs wiener estimation for reflectance from images and returns an sRGB reconstucted image
%Input arguments
%g = the image matrix g(filter band, x , y , image rgb values)
%spectrum = the measured originalReflectance spectrum of the sample
%id = the point id information
%
%Optional input arguments
%'pixelValueSelectionMethod' = {'green', 'rms', 'adjusted', 'extended', 'unchanged'} for selecting a single value from the 3 RGB values
%'smoothingMatrixMethod' = {'markovian 0.99', 'KCor all data', 'KCor macbeth', 'KCor all specimen', ...
%                'KCor all same malignancy', 'KCor all same fixation' 'KCor same malignancy', 'KCor same malignancy, fixation'}; for selecting the smoothing matrix
%'rho' =  parameter for markovian smoothing matrix , default: 0.99
%'variance' =  parameter for noise covariance matrix , default: 1
%'noiseType' = to include noise component,{'none', 'independent', 'spatial', 'givenSNR', 'white gaussian', 'fromOlympus'}
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
    mask = g.Mask;
    g = g.Patch; 
else 
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

if isfield(options, 'variance')
    variance = options.variance;
else 
    variance = m.noiseparam';
end

if strcmp(options.pixelValueSelectionMethod, 'rgb')
    coeff = ones(3, 1) / 10^3;
else
    c = matfile(fullfile(options.systemdir, 'coeff.mat'));
    pixelValueSelectionMethods = {'green', 'rms', 'adjusted', 'extended', 'rgb'};
    coeff = squeeze(c.coeff(id.Index, min(find(strcmp(pixelValueSelectionMethods, options.pixelValueSelectionMethod)), 3), :))';
    clear('c');
end

if strcmp(pixelValueSelectionMethod, 'extended')
    if ~isempty(coeff)
        coeff = [coeff(1:2), coeff(2), coeff(3:5), coeff(6), coeff(6:7)];
    end
    variance = [variance(1:2), variance(2), variance(3:5), variance(6), variance(6:7)];
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

G = valueSelect(g, pixelValueSelectionMethod); % convert from 4D MSI+rgb to 3D MSI+grey
[msibands, height, width] = size(G);

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
        M = m.Cor_Sample(id.SampleId);
        
    case 'Cor_Malignancy'
        if isBenign
            M = m.Cor_Benign;
        else
            M = m.Cor_Malignant;
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
            M = m.Cor_SampleBenign(id.SampleId);
        else
            M = m.Cor_SampleMalignant(id.SampleId);
        end
        
    case 'Cor_SampleMalignancyFixing' % average xcorr of all fixed cancer measured data for fixed cancer sample
        if isBenign && strcmp(fixing, 'Cut')
            M = m.Cor_SampleBenignCut(id.SampleId);
        elseif isBenign && strcmp(fixing, 'Fixed')
            M = m.Cor_SampleBenignFixed(id.SampleId);
        elseif isBenign && strcmp(fixing, 'Unfixed')
            M = m.Cor_SampleBenignUnfixed(id.SampleId);
        elseif ~isBenign && strcmp(fixing, 'Cut')
            M = m.Cor_SampleMalignantCut(id.SampleId);
        elseif ~isBenign && strcmp(fixing, 'Fixed')
            M = m.Cor_SampleMalignantFixed(id.SampleId);
        elseif ~isBenign && strcmp(fixing, 'Unfixed')
            M = m.Cor_SampleMalignantUnfixed(id.SampleId);
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
if contains(noiseType, 'independent')
    attr = strsplit(noiseType, {' ', '^{', '}'});
    if (numel(attr) > 1)
        options.noiseOrder = str2double(attr{2})^str2double(attr{3});
    end
    Kn = diag(options.noiseOrder*ones(size(variance)));

elseif contains(noiseType, 'givenSNR')
    defaultSNR = 36; %dB
    attr = strsplit(noiseType, {' ', 'dB'});
    if (numel(attr) > 1)
        snr = str2double(attr{2});
    elseif isfield(options, 'snr')
        snr = options.snr;
    else
        snr = defaultSNR; %dB
    end
    variance = (trace(HMH) / (msibands * 10^(snr / 10))) * ones(msibands, 1);
    Kn = diag(variance); %different at every channel def. 0.0379
    
elseif contains(noiseType, 'white gaussian')
    attr = strsplit(options.noiseType, {' ', '^{', '}'});
    if (numel(attr) > 1)
        options.noiseOrder = str2double(attr{3})^str2double(attr{4});
    end
    Kn = diag( randn(1, msibands) .* options.noiseOrder); %different at every channel

elseif strcmp(noiseType, 'fromOlympus')
    Kn = diag(variance); %different at every channel
    
elseif strcmp(noiseType, 'none')
    Kn = 0;
    
elseif strcmp(noiseType, 'spatial') || contains(noiseType, 'spatial')
    if isfield(options, 'windowDim')
        windowDim = options.windowDim;
    else 
        windowDim = 3;
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
    
    if exist('sigma1', 'var') && ~exist('sigma2', 'var')
        variance = ones(msibands,1) * (sqrt(0.5) * sigma1 + sigma2)^2; %else use noise variance from Olympus 
    end
    
    Kn = diag(variance);
    
else 
    error('Not implemented');
end
    
% Covariance matrix of the additive noise ends

Gres = reshape(G, msibands, height*width); % msi = im2double(diag(coeff) * double(im2uint16(G(:,row,col))) );
mask = reshape(mask, 1, height*width);
    
if strcmp(noiseType, 'spatial')
%% Perform Spatially Adaptive Wiener estimation for all pixels in an image area

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
    estimatedReflectance = reshape(estimatedReflectance, numel(wavelength), height, width);
else   
    
    if ~isempty(spectrum)    
        [rmse, rmseIdx] = min(Rmse(spectrum, estimatedReflectance));
        [nmse, nmseIdx] = min(Nmse(spectrum, estimatedReflectance));
        if (rmseIdx ~= nmseIdx)
            disp('Different minimum location.')
        end
        estimatedReflectance = estimatedReflectance(:,rmseIdx);
    end
end

% Perform the estimation for all pixels in an image area ends

end

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