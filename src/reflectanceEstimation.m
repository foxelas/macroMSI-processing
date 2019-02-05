function [estimatedReflectance, rmse, nmse, minIdx] = reflectanceEstimation(MSI, mask, spectrum, id, options)

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

if isempty(mask)
     mask = ones(size(MSI, 2), size(MSI, 3));
end

if size(spectrum,2) ~= 1 
    spectrum = spectrum';
end

precomputedFile = fullfile(options.systemdir, 'precomputedParams.mat'); %pre-set parameters
smoothingMatrixMethod = options.smoothingMatrixMethod;
pixelValueSelectionMethod = options.pixelValueSelectionMethod;
noiseType = options.noiseType;

if strcmp(smoothingMatrixMethod, 'markovian')
    defaultRho = 0.99;
    if isfield(options, 'rho')
        rho = options.rho;
    else 
        rho = defaultRho;
    end
end

if strcmp(smoothingMatrixMethod, 'adaptive')
    defaultAlpha = 0.7;
    if isfield(options, 'alpha')
        alpha = options.alpha;
    else
        alpha = defaultAlpha;
    end
    
    defaultGamma = 2;
    if isfield(options, 'gamma')
        gamma = options.gamma;
    else
        gamma = defaultGamma;
    end
end 

signalCov = [0.0290584774952317;0.0355809665018991;0.0254338785693633;0.0278002826809506;0.00535859811159383;0.0263030884980115;0.0198921193827871];

if strcmp(options.pixelValueSelectionMethod, 'rgb')
    coeff = ones(1,3) / 10^3;
else
    load(precomputedFile, 'coeff');
    pixelValueSelectionMethods = {'green', 'rms', 'adjusted', 'extended', 'rgb'};
    pvmIdx = min(find(strcmp(pixelValueSelectionMethods, options.pixelValueSelectionMethod)), 3);
    coeff = squeeze(Coefficients(id.Index, pvmIdx, 1:7))'; %id.CoeffIndex
end


if strcmp(pixelValueSelectionMethod, 'extended')
    load(precomputedFile, 'Hextended');
    H = Hextended;
    if ~isempty(coeff)
        coeff = [coeff(1:2), coeff(2), coeff(3:5), coeff(6), coeff(6:7)];
    end
elseif strcmp(pixelValueSelectionMethod, 'green')
    load(precomputedFile, 'Hgreen');
    H = Hgreen;
    
elseif strcmp(pixelValueSelectionMethod, 'rms')
    load(precomputedFile, 'Hrms');
    H = Hrms;

elseif strcmp(pixelValueSelectionMethod, 'adjusted') 
    load(precomputedFile, 'Hadjusted');
    H = Hadjusted; 
    
elseif strcmp(pixelValueSelectionMethod, 'rgb') 
    load(precomputedFile, 'Hrgb');
    H = Hrgb; 
else
    error('Unsupported "pixelValueSelectionMethod".');
end

coeff = coeff * 5; % to adapt coefficients from bandwidth 1nm to bandwidth 5nm (Source is 1x141, esimation is 1x81)
H = diag(coeff) * H; % illumination x sensitivity

if isfield(options, 'SVDTol')
    hasSVDTol = true;
    tol = options.SVDTol;
else
    hasSVDTol = false;
    tol = 1;
end

% Argument parsing ends

%% Generate smoothing matrix for estimation

G = raw2msi(MSI, pixelValueSelectionMethod); % convert from 4D MSI+rgb to 3D MSI+grey
[msibands, height, width] = size(G);
%G = G ./ options.luminanceCorrection;

% computed beforehand with prepareSmoothingMatrix
isBenign = id.IsNormal;
if (id.IsCut)
    fixing = 'Cut';
elseif (id.IsFixed)
    fixing = 'Fixed';
else
    fixing = 'Unfixed';
end

switch smoothingMatrixMethod
    case {'Markovian 0.99', 'markovian'}
        if rho ~= defaultRho
            error('Invalid correlation parameter for the markovian process');
        end      
        load(precomputedFile, 'Cor_Markovian');
        M = mean(signalCov) * Cor_Markovian;
        
    case 'Cor_All' % average xcorr of all measured data
        load(precomputedFile, 'Cor_All');
        M = Cor_All; 
        
    case 'Cor_Macbeth' % average xcorr of all macbeth data
        load(precomputedFile, 'Cor_Macbeth');
        M = Cor_Macbeth; 
        
    case 'Cor_Sample' % average xcorr of all measured data for each sample
        load(precomputedFile, 'Cor_Sample');
        M = squeeze(Cor_Sample(id.SampleId,:,:));
        
    case 'Cor_Malignancy'
        if isBenign
            load(precomputedFile, 'Cor_Benign');
            M = Cor_Benign; %m.Mn; 
        else
            load(precomputedFile, 'Cor_Malignant');
            M = Cor_Malignant; %m.Mc;
        end
        
    case 'Cor_Fixing'
        if strcmp(fixing, 'Cut')
            load(precomputedFile, 'Cor_Cut');
            M = Cor_Cut;
        elseif strcmp(fixing, 'Fixed')
            load(precomputedFile, 'Cor_Fixed');
            M = Cor_Fixed;
        else
            load(precomputedFile, 'Cor_Unfixed');
            M = Cor_Unfixed;
        end
        
    case 'Cor_SampleMalignancy' % average xcorr of all measured cancerous data for cancerous sample
        if isBenign
            load(precomputedFile, 'Cor_SampleBenign');
            M = squeeze(Cor_SampleBenign(id.SampleId,:,:));
        else
            load(precomputedFile, 'Cor_SampleMalignant');
            M = squeeze(Cor_SampleMalignant(id.SampleId,:,:));
        end
        
    case 'Cor_SampleMalignancyFixing' % average xcorr of all fixed cancer measured data for fixed cancer sample
        if isBenign && strcmp(fixing, 'Cut')
            load(precomputedFile, 'Cor_SampleBenignCut');
            M = squeeze(Cor_SampleBenignCut(id.SampleId,:,:));
        elseif isBenign && strcmp(fixing, 'Fixed')
            load(precomputedFile, 'Cor_SampleBenignFixed');
            M = squeeze(Cor_SampleBenignFixed(id.SampleId,:,:));
        elseif isBenign && strcmp(fixing, 'Unfixed')
            load(precomputedFile, 'Cor_SampleBenignUnfixed');
            M = squeeze(Cor_SampleBenignUnfixed(id.SampleId,:,:));
        elseif ~isBenign && strcmp(fixing, 'Cut')
            load(precomputedFile, 'Cor_SampleMalignantCut');
            M = squeeze(Cor_SampleMalignantCut(id.SampleId,:,:));
        elseif ~isBenign && strcmp(fixing, 'Fixed')
            load(precomputedFile, 'Cor_SampleMalignantFixed');
            M = squeeze(Cor_SampleMalignantFixed(id.SampleId,:,:));
        elseif ~isBenign && strcmp(fixing, 'Unfixed')
            load(precomputedFile, 'Cor_SampleMalignantUnfixed');
            M = squeeze(Cor_SampleMalignantUnfixed(id.SampleId,:,:));
        else 
            error('Unsupported type')
        end
        
    case 'Cor_MalignancyFixing'
        if isBenign && strcmp(fixing, 'Cut')
            load(precomputedFile, 'Cor_BenignCut');
            M = Cor_BenignCut;
        elseif isBenign && strcmp(fixing, 'Fixed')
            load(precomputedFile, 'Cor_BenignFixed');
            M = Cor_BenignFixed;
        elseif isBenign && strcmp(fixing, 'Unfixed')
            load(precomputedFile, 'Cor_BenignUnfixed');
            M = Cor_BenignUnfixed;
        elseif ~isBenign && strcmp(fixing, 'Cut')
            load(precomputedFile, 'Cor_MalignantCut');
            M = Cor_MalignantCut;
        elseif ~isBenign && strcmp(fixing, 'Fixed')
            load(precomputedFile, 'Cor_MalignantFixed');
            M = Cor_MalignantFixed;
        elseif ~isBenign && strcmp(fixing, 'Unfixed')
            load(precomputedFile, 'Cor_MalignantUnfixed');
            M = Cor_MalignantUnfixed;
        else 
            error('Unsupported type')
        end
        
    case 'adaptive'
        %% Based on "Reflectance reconstruction for multispectral imaging by adaptive Wiener estimation"[Shen2007]
        options.smoothingMatrixMethod = 'Cor_Sample';
        rhat = reflectanceEstimation(MSI, mask, spectrum, id, options);
        M = adaptiveSmoothingMatrix(rhat, options.systemdir, alpha, gamma);
       
    otherwise
        error('Unexpected smoothing matrix method. Abort execution.')
end
% Generate smoothing matrix for estimation ends

%% Wiener Matrices 
MH = M * H';
HMH = H * M * H';

%% Covariance matrix of the additive noise
noiseParts = strsplit(noiseType, {' ', ','});
if (numel(noiseParts) > 1); options.noiseParam = cellfun(@str2double , noiseParts(2:end)); end
hasNoiseParam = isfield(options, 'noiseParam');

if contains(noiseType, 'sameForChannel')
    if (hasNoiseParam); variance = options.noiseParam * ones(msibands, 1); else; variance = 0.001 * ones(msibands, 1); end
    
elseif contains(noiseType, 'diffForChannel')
    if (hasNoiseParam); variance = options.noiseParam; else; variance = [0.0031, 0.0033, 0.0030, 0.0031, 0.0032, 0.0029, 0.0024]; end
    
elseif contains(noiseType, 'givenSNR')
    if (hasNoiseParam); variance = (trace(HMH) / (msibands * 10^(options.noiseParam / 10))) * ones(msibands, 1); else; variance = (trace(HMH) / (msibands * 10^( 17 / 10))) * ones(msibands, 1); end

elseif contains(noiseType, 'white gaussian')
    if (hasNoiseParam); variance = (randn(1, msibands) .* options.noiseParam).^2; else; variance = (randn(1, msibands) .* 0.0001).^2; end

elseif strcmp(noiseType, 'fromOlympus')
    noiseparam = [1.62215000000000e-05;1.57000000000000e-05;7.55000000000000e-06;5.03000000000000e-06;8.38000000000000e-06;0.000148000000000000;1.48000000000000e-05];
    variance =  (noiseparam').^2 ;
    
elseif strcmp(noiseType, 'none')
    variance = zeros(1, msibands);
    
elseif strcmp(noiseType, 'spatial') || contains(noiseType, 'spatial')
    if isfield(options, 'windowDim'); windowDim = options.windowDim; else; windowDim = 5; end
    [windowKernel, windowElements] = makeKernel(windowDim, height, width);  
    if (hasNoiseParam); variance = ones(msibands,1) * (sqrt(0.5) * options.noiseParam(1) + options.noiseParam(2))^2; else; variance = ones(msibands,1) * (sqrt(0.5) * 0.001 + 0.03)^2; end
        
else 
    error('Not implemented');
end

if (length(variance) < msibands)
    variance = [variance(1:2), variance(2), variance(3:5), variance(6), variance(6:7)];
elseif (length(variance) > msibands)
    variance = [variance(1), variance(4), variance(7)];    
end
Kn = diag(variance); %different at every channel
% Covariance matrix of the additive noise ends

hw = height * width;
Gres = reshape(G, msibands, hw); % msi = im2double(diag(coeff) * double(im2uint16(G(:,row,col))) );
activeRegionIdx = sub2ind([height,width], find(mask));
   
if strcmp(noiseType, 'spatial')
%% Perform Spatially Adaptive Wiener estimation for all pixels in an image area
%  Based on "A Spatially AdaptiveWiener Filter for Reflectance Estimation"[Urban2008]

    meanG = zeros(msibands, height, width);
    centeredG = zeros(msibands, height, width);
    Kcov = zeros(msibands, height, width);
    for b = 1:msibands
        Gb = squeeze(G(b,:,:));
        meanG(b,:,:) = conv2(Gb, windowKernel ./ windowElements,'same');        
        centeredG(b,:,:) = reshape(Gb, [1, height, width]) - meanG(b,:,:);
        sdGb = stdfilt(Gb, windowKernel);
        Kcov(b, :, :) = sdGb.^2;
    end
    
    estimatedReflectance = zeros(length(spectrum), length(activeRegionIdx));
    for p = 1:length(activeRegionIdx)
        [i, j] = ind2sub([height,width], activeRegionIdx(p));
        meanGij = meanG(:,i,j);
        centeredGij = centeredG(:,i,j);
        Kij = diag(Kcov(:,i,j));
        W = multiplyToInverse(Kij, Kij + Kn); % same as kInWindow * inv(kInWindow + Kn)
        What = multiplyToInverse(MH, HMH + Kij - W * Kij, hasSVDTol, tol);
        estimatedReflectance(:, p) = What * W * centeredGij + What * meanGij;  
    end
    
else    
    div = multiplyToInverse(MH , HMH + Kn, hasSVDTol, tol); % M 401x401, H' 401x7, inv() 7x7
    estimatedReflectance = div * Gres(:, activeRegionIdx); % 401 x 100 
end  

%% Show all estimates in the region
if (false)
    figure(3);
    clf(3);
    [rmse, rmseIdx] = min(Rmse(spectrum, estimatedReflectance));
    bb = estimatedReflectance(:,rmseIdx);
    [mini, minj] = ind2sub([height, width], activeRegionIdx(rmseIdx));
    minIdx = [mini, minj];
    hold on
    for kkk = 1:size(estimatedReflectance,2)
    plot(estimatedReflectance(:,kkk));
    end
    plot(spectrum, 'm*')
    plot(bb, 'g*')
    scatter( ([450, 465, 505, 525, 575, 605, 630] - 380) /5, squeeze( raw2msi(MSI(:, minIdx(1), minIdx(2), :), 'adjusted')) , 'bo');
    hold off
    pause(0.2)
end

%% Perform Wiener estimation for all pixels in an image area
if (height > 200 ||  width > 200) %in this case the mask is ones(height, width)
    estimatedReflectance = max(estimatedReflectance, 0);
    estimatedReflectance = min(estimatedReflectance, 1);
    
    estimatedReflectance = reshape(estimatedReflectance, numel(spectrum), height, width);
else   

    idx = any(estimatedReflectance < 0) | any(estimatedReflectance > 1);
    if any(estimatedReflectance(:) < 0 ) || any(estimatedReflectance(:) > 1 )
        estimatedReflectance = max(estimatedReflectance, 0);
        estimatedReflectance = min(estimatedReflectance, 1);
        warning('Estimated reflectance spectrum is out of bounds for all pixels in the region!(ID index = %d)', id.Index);
    else
        estimatedReflectance = estimatedReflectance(:,~idx);
    end
    
    [rmse, rmseIdx] = min(Rmse(spectrum, estimatedReflectance));
    [nmse, ~] = min(Nmse(spectrum, estimatedReflectance));
    estimatedReflectance = estimatedReflectance(:,rmseIdx);
    [mini, minj] = ind2sub([height, width], activeRegionIdx(rmseIdx));
    minIdx = [mini, minj]; 
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

% function di = reflectanceDistance(ri, rhat, a)
%     di = a * mean(abs(ri / norm(ri) - rhat / norm(rhat)))...
%         + (1-a) * max(abs(ri / norm(ri) - rhat / norm(rhat)));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function k = replicationTimes(d, dmax, gamma)
    if (nargin < 3)
        gamma = 1;
    end
    k = floor((dmax / d)^gamma + 0.5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function adaptedM = adaptiveSmoothingMatrix(rhat, systemdir, a, gamma)

    load(fullfile(systemdir, 'in.mat'), 'Spectra', 'SpectraNames');
    [~, idxs] = unique(SpectraNames);
    r = num2cell( Spectra(idxs,:));
    d = cellfun(@(x) DiscreteFrechetDist(x, rhat), r); % or reflectanceDistance
    reps = arrayfun(@(x) replicationTimes(x, max(d), gamma), d);
    spectra = zeros(length(rhat), sum(reps));
    j = 0;
    for i = 1:length(d)
        k = reps(i);
        spectra(:, (j+1):(j+k)) = repmat(r{i}, 1, k);
        j = j + k;
    end

    means = mean(spectra,2);
    adaptedM = 1 / (size(spectra, 2) - 1) * (spectra - means) * (spectra - means)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [windowKernel, windowElements] = makeKernel(windowDim, height, width)

    windowKernel = ones(windowDim);
    windowElements = windowDim^2;
    
    if size(windowKernel, 1) > height || size(windowKernel,2) > width
        mindim = min(height, width);
        if mod(mindim,2)==0
            mindim = mindim - 1;
        end
        windowKernel = ones(mindim);
        windowElements = min(height, width);
    end
end