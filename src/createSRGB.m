function [sRGB, Lab16] = createSRGB(I, method, outName, id, options, adaptationModel)
%CREATESGB generates the sRGB based on the MSI reflectances
%   Use: sRGB = createSRGB( g, 'out.jpg', [50, 465, 505, 525, 575, 605,
%   630], 'd65');
%   reflectance: the input matrix with size MxNx7, if it's 4D then converts
%   to MxNx7 using 'adjusted'
%   outName: if it exists, then the sRGB image is saved
%   adaptationModel: creates the conversion matrixn from XYZ to sRGB based
%   on an color adaptation model

if strcmp(method, 'original')
    
    wavelength = 380:5:780;
    [w, r, c, ~] = size(I);
    
    % load color matching functions
    [lambdaMatch, xFcn, yFcn, zFcn] = colorMatchFcn('1931_FULL');
    CMF = interp1(lambdaMatch, [xFcn; yFcn; zFcn]', wavelength', 'pchip', 0);
    
    % Illumination
    [lambdaMatch, d65] = illuminant('d65');
    d65 = interp1(lambdaMatch, d65, wavelength', 'nearest');
    
    XYZ = bsxfun(@times, CMF, d65);
    
    %Compute the normalisation factor k
    k = 100 ./ sum(XYZ(:, 2));
    
    options.smoothingMatrixMethod = 'KCor all specimen';
    options.pixelValueSelectionMethod = 'extended';
    options.noiseType = 'givenSNR';
    options.snr = 20;
    
    reflectance = reflectanceEstimation(I, [], id, options);
    reflectance = permute(reflectance, [2, 3, 1]);
    reflectance = reshape(reflectance, r*c, numel(wavelength));
    XYZrecon = reflectance * XYZ;
    XYZrecon = reshape(XYZrecon, r, c, 3);
    XYZrecon = k * XYZrecon ./ 100;
    
    % Convert the XYZ values to sRGB.
    
    C = makecform('XYZ2Lab');
    Lab = applycform(XYZrecon, C);
    Lab16 = uint16(Lab./100.*((2^16) - 1));
    
    C = makecform('XYZ2srgb');
    sRGB = applycform(XYZrecon, C);
    imwrite(sRGB, strcat(options.systemdir, 'sRGB.tif'))
    
    
elseif strcmp(method, 'tutorial')
    g = valueSelect(I, 'adjusted');
    [w, r, c] = size(g);
    if (w == 7)
        wavelength = [450, 465, 505, 525, 575, 605, 630]; % 7 BANDS
    elseif (w == 3)
        wavelength = [465, 535, 595]; % based on the sensitivity values
    end
    
    % load color matching functions
    [lambdaMatch, xFcn, yFcn, zFcn] = colorMatchFcn('1931_FULL');
    
    % Interpolate the input wavelength in the color matching functions. 7x3
    xyzbar = interp1(lambdaMatch', [xFcn; yFcn; zFcn]', wavelength', 'pchip', 0);
    
    % Illumination
    [lambda, d65] = illuminant('d65');
    d65 = interp1(lambda, d65, wavelength', 'nearest');
    
    % Radiance
    radiances = bsxfun(@times, reshape(g, [r * c, w]), d65');
    
    XYZ = (xyzbar' * radiances')';
    XYZ = reshape(XYZ, [r * c, 3]);
    XYZ = max(XYZ, 0);
    XYZ = XYZ / max(XYZ(:)); %normalization [0,1]
    
    if ~isempty(adaptationModel)
        targetXYZ = [0.95047; 1.00000; 1.08883]; %d65 white
        sourceXYZ = [0.0725; 0.0418; 0.1583];
        M = cbCAT(sourceXYZ, targetXYZ, adaptationModel);
        
    else
        % convert this XYZ representation of the image to the default RGB colour representation sRGB (IEC_61966-2-1),
        % Forward transformation from 1931 CIE XYZ values to sRGB values (Eqn 6 in
        % IEC_61966-2-1.pdf).
        M = [3.2406, -1.5372, -0.4986; ...
            -0.9689, 1.8758, 0.0414; ...
            0.0557, -0.2040, 1.0570];
    end
    
    sRGB = (M * XYZ')';
    
    % Gamma correction
    a = 0.055;
    t = 0.0031308;
    gamma = 2.4;
    
    for i = 1:r * c
        for j = 1:3
            if (sRGB(i, j) <= t)
                sRGB(i, j) = 12.92 * sRGB(i, j);
            else
                sRGB(i, j) = (1 + a) * sRGB(i, j)^(1 / gamma) - a;
            end
        end
    end
    
    % Reshape to recover shape of original input.
    sRGB = reshape(sRGB, [r, c, 3]);
    
    %clip values to [0,1]
    sRGB = max(sRGB, 0);
    sRGB = min(sRGB, 1);
    
    figure(1);
    imshow(sRGB, 'Border', 'tight');
    title('sRGB image resulting from the MSI image')
    
    % Plot Chromacity chart
    x = bsxfun(@(x, s) x./s, XYZ(:, 1), sum(XYZ, 2));
    y = bsxfun(@(x, s) x./s, XYZ(:, 2), sum(XYZ, 2));
    figure(2);
    cieplot();
    hold on
    scatter(x, y, 10);
    hold off
else
    warning('Not implemented')
end

figure(1);
imshow(sRGB);
title('sRGB from MSI data')

if ~(isempty(outName))
    imwrite(sRGB, strcat(outName, '_sRGB', '.tif'), 'tif');
    
    %Write to TIFF (in the format that Photoshop can read the L,a,b channels)
    imwrite(Lab16(:, :, 1), strcat(outName, '_L', '.tif'), 'tif');
    imwrite(Lab16(:, :, 2)./2+((2^16) / 2), strcat(outName, '_a', '.tif'), 'tif');
    imwrite(Lab16(:, :, 3)./2+((2^16) / 2), strcat(outName, '_b', '.tif'), 'tif');
end

%% Makes a 2D histogram style image of wavelength versus intensity of the spectral map

I = valueSelect(I, 'adjusted');

vec = zeros(w, r*c);
vech = zeros(256, w);
for i = 1:w
    vec(i, :) = reshape(squeeze(I(i, :, :)), 1, r*c);
    h = hist(vec(i, :).*255, 0:1:255);
    vech(1:256, i) = h;
end
vech = flipud(vech*(2^16 - 1)/max(max(vech)));
v = figure(2);
imagesc(vech);
xlabel('Frequency band')
ylabel('Image intensity value')
title('Intensity Histogram')
colorbar

if ~(isempty(outName))
    saveas(v, strcat(outName, '_binnedspectrum.jpg'));
end

end
