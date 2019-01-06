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
%'windowDim' = the dimension of square pixel neighborhood (odd number), priority: high, default: 3
%'windowDiam' = the diameter of round pixel neighborhood (odd number), priority: medium
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

m = matfile(fullfile(options.systemdir, 'param.mat')); %pre-set parameters

wavelengthN = size(sensitivity, 1);
wavelength = linspace(380, 780, wavelengthN);

if isfield(options, 'smoothingMatrixMethod')
    smoothingMatrixMethod = options.smoothingMatrixMethod;
else 
    smoothingMatrixMethod = 'markovian';
end

if strcmp(smoothingMatrixMethod, 'markovian')
    if isfield(options, 'rho')
        rho = options.rho;
    else 
        rho = 0.99;
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

if isfield(options, 'isSpatial')
    isSpatial = options.isSpatial;
    if isSpatial 
        if isfield(options, 'windowDim')
            windowDim = options.windowDim;
            windowType = 'square';
            windowKernel = ones(windowDim);
            windowElements = (2 * windowDim + 1)^2;

        elseif isfield(options, 'windowDiam')
            windowDiam = options.windowDiam;
            windowType = 'round';
            radius = ceil(windowDiam / 2);
            centerY = radius;
            centerX = radius;
            [cols, rows] = meshgrid(1:windowDiam, 1:windowDiam);
            windowKernel = (rows - centerY).^2 + (cols - centerX).^2 <= radius.^2;
            windowElements = sum(windowKernel(:));
            clear('centerX', 'centerY', 'radius', 'cols', 'rows');
        else 
            windowDim = 3;
            windowType = 'square';
            windowKernel = ones(windowDim);
            windowElements = (2 * windowDim + 1)^2;
        end
    end
else 
    isSpatial = false
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

if ~isempty(id)
    sample = generateName([], 'sample', [], id);
    
    if (id.IsNormal)
        cancerType = 'n';
    else
        cancerType = 'c';
    end
    
    if (id.IsCut)
        fixType = 'c';
    elseif (id.IsFixed)
        fixType = 'f';
    else
        fixType = 'u';
    end
end

% Argument parsing ends

%% Generate smoothing matrix for estimation

G = valueSelect(g, pixelValueSelectionMethod); % convert from 4D MSI+rgb to 3D MSI+grey
[msibands, height, width] = size(G);

coeff = coeff * 5; % to adapt coefficients from bandwidth 1nm to bandwidth 5nm (Source is 1x141, esimation is 1x81)

% computed beforehand with prepareSmoothingMatrix
switch smoothingMatrixMethod
    case {'Markovian 0.99', 'markovian'}
        M = m.M_markov_99;
        M = mean(m.signalCov) * M;
        
    case 'corr all spectra' % average xcorr of all measured data
        M = m.M_measured_xcorr';
        
    case 'corr macbeth spectra' % average xcorr of all macbeth data
        M = m.M_measured_xcorr_macbeth;
        
    case 'corr sample spectra' % average xcorr of all measured data for each sample
        M = eval(strcat('m.M_', sample));
        
    case 'corr same malignancy all spectra'
        M = eval(strcat('m.M', cancerType));
        
    case 'corr same fixing all spectra'
        M = eval(strcat('m.M', '_', fixType));
        
    case 'corr same malignancy sample spectra' % average xcorr of all measured cancerous data for cancerous sample
        M = eval(strcat('m.M', cancerType, '_', sample));
        
    case 'KCor same malignancy, fixation' % average xcorr of all fixed cancer measured data for fixed cancer sample
        M = eval(strcat('m.M', cancerType, fixType, '_', sample));        
        
    case 'spatiospectral'
        %unfinished
        spectralcor = eval(strcat('m.M', cancerType, '_', sample));
        
        rho = 0.985;
        rho1 = rho;
        rho2 = rho;
        R1 = markovianMatrix(rho1);
        R2 = markovianMatrix(rho2);
        spatialcor = sigma * sigma * kron(R1, R2);
        
        
    case 'obi2001' %%unfinished
        S = valueSelect(sensitivity, pixelValueSelectionMethod, wavelength);
        D = S * illumination * repmat(spectrum', 7, 1);
        A = spectrum' * spectrum;
        Aplus = (A' * A) \ (A'); % A+ = (A'*A)^-1*A'   401x401
        B = spectrum * spectrum'; % should be changed to principal component  401xR
        P = Aplus * B; % Ks = D * P
        P0 = eye(wavelengthN);
        %rms =@(x) (norm(P - x)^ 2); %function similar to yours
        P = fminsearch(rms, P0); %find the x
        M = D * P;
        
    otherwise
        error('Unexpected smoothing matrix method. Abort execution.')
end

% Generate smoothing matrix for estimation ends

%% Illumination X sensitivity

H = illumXsensitivity(illumination, sensitivity, pixelValueSelectionMethod); % illumination * sensitivity
H = diag(coeff) * H;

% Illumination X sensitivity ends

%% Covariance matrix of the additive noise
if contains(noiseType, 'independent')
    attr = strsplit(noiseType, {' ', '^{', '}'});
    if (numel(attr) > 1)
        options.noiseOrder = str2double(attr{2})^str2double(attr{3});
    end
    Kn = diag(options.noiseOrder*ones(size(variance)));

elseif contains(noiseType, 'givenSNR')
    attr = strsplit(noiseType, {' ', 'dB'});
    if (numel(attr) > 1)
        snr = str2double(attr{2});
    elseif isfield(options, 'snr')
        snr = options.snr;
    else
        snr = m.snr; %dB
    end
    variance = (trace(H*M*H') / (msibands * 10^(snr / 10))) * ones(msibands, 1);
    Kn = diag(variance); %different at every channel
    
elseif contains(noiseType, 'white gaussian')
    attr = strsplit(options.noiseType, {' ', '^{', '}'});
    if (numel(attr) > 1)
        options.noiseOrder = str2double(attr{3})^str2double(attr{4});
    end
    Kn = diag(randn(1, length(coeff)).*options.noiseOrder); %different at every channel

elseif strcmp(noiseType, 'fromOlympus')
    Kn = diag(variance); %different at every channel
    
elseif strcmp(noiseType, 'spatial')
    if isfield(noiseType, 'sigma1')
        sigma1 = options.sigma1;
    else
        sigma1 = 0.002; %0.002 to 0.01
    end
    if isfield(noiseType, 'sigma2')
        sigma2 = options.sigma2;
    else
        sigma2 = 0.002; %0.002 to 0.01
    end
    
    variance = ones(msibands,1) * (sqrt(0.5) * sigma1 + sigma2)^2;
    Kn = diag(variance);
    
elseif strcmp(noiseType, 'none')
    Kn = 0;
    
elseif strcmp(noiseType, 'dependent')
    error('not implemented.')
    
else 
    error('Not implemented');
end
    
% Covariance matrix of the additive noise ends

%% Wiener matrix
div = M * H' / (H * M * H' + Kn); % M 401x401, H' 401x7, inv() 7x7

if isfield(options, 'withSVDinHMHKn')
    tol = 1; % singular value threshold
    div = M * H' * pinv(H*M*H'+Kn, tol);
end
%Wiener matrix ends

if isSpatial
%% Perform Spatially Adaptive Wiener estimation for all pixels in an image area
    meansInWindow = conv2(G, windowKernel ./ windowElements,'same');
    kInWindow = conv2((G - meansInWindow)*(G - meansInWindow)', windowKernel ./ (windowElements - 1), 'same');
    W = kInWindow / (kInWindow + Kn); % same as kInWindow * inv(kInWindow + Kn)
    What = M * H' / (H * M * H' + kInWindow - W * kInWindow);
    div = What * W * (G - meansInWindow) + What * meansInWindow;
    div = reshape(div, msibands, height*width); %             msi = im2double(diag(coeff) * double(im2uint16(G(:,row,col))) );
end  

%% Perform Wiener estimation for all pixels in an image area
G = reshape(G, msibands, height*width); %             msi = im2double(diag(coeff) * double(im2uint16(G(:,row,col))) );
estimatedReflectance = div * G; % 401 x 100


if (size(g, 2) > 30)
    estimatedReflectance = reshape(estimatedReflectance, numel(wavelength), height, width);
else
    estimatedReflectance = mean(estimatedReflectance, 2);
    
    if ~isempty(spectrum)
        r = spectrum;
        rhat = estimatedReflectance;
        N = length(r);
        
        % Root Mean Square Error
        rmse = sqrt((r - rhat)'*(r - rhat)/N);

        % Normalized Mean Square Error
        nmse = (r - rhat)' * (r - rhat) / (r' * r);
    end
end

% Perform the estimation for all pixels in an image area ends

end
