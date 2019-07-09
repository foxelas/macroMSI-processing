function [sRGB, Lab16] = createSRGB(I, method, id, options, adaptationModel)
%CREATESGB generates the sRGB based on the MSI reflectances
%   Use: sRGB = createSRGB( g, 'out.jpg', [50, 465, 505, 525, 575, 605,
%   630], 'd65');
%   reflectance: the input matrix with size MxNx7, if it's 4D then converts
%   to MxNx7 using 'adjusted'
%   outName: if it exists, then the sRGB image is saved
%   adaptationModel: creates the conversion matrixn from XYZ to sRGB based
%   on an color adaptation model

wavelength = 380:5:780;
[w, r, c, ~] = size(I);


% load color matching functions
[lambdaMatch, xFcn, yFcn, zFcn] = colorMatchFcn('1931_FULL');
CMF = interp1(lambdaMatch, [xFcn; yFcn; zFcn]', wavelength', 'pchip', 0);

% Illumination
[lambdaMatch, d65] = illuminant('d65');
d65 = interp1(lambdaMatch, d65, wavelength', 'nearest');

options.smoothingMatrixMethod = 'Cor_All';
options.pixelValueSelectionMethod = 'extended';            
options.noiseType = 'spatiospectralolympus';
options.rho = 0.6;
options.windowDim = 3;
options.noiseParam = 0.0001;   

reflectance = reflectanceEstimation(I, [], [], id, options);
reflectance = permute(reflectance, [2, 3, 1]);

if strcmp(method, 'original')
       
    XYZ = bsxfun(@times, CMF, d65);
    
    %Compute the normalisation factor k
    reflectance = reshape(reflectance, r*c, numel(wavelength));
    XYZrecon = reflectance * XYZ;
    XYZrecon = reshape(XYZrecon, [r, c, 3]);
    %k = 1 / sum(XYZ(:, 2));
    %XYZrecon = k * XYZrecon;
    XYZrecon =  XYZrecon / max(XYZrecon(:));
    
    % Convert the XYZ values to sRGB.
    
    C = makecform('XYZ2Lab');
    Lab = applycform(XYZrecon, C);
    Lab16 = uint16(Lab./100.*((2^16) - 1));
    
    sRGB = xyz2rgb(XYZrecon, 'WhitePoint','d65');

    sRGB2 = reshape(xyz2srgbCustom( reshape(XYZrecon , r*c, 3)), [r,c,3]);
    
elseif strcmp(method, 'tutorial')
    
    % Radiance
    radiances = bsxfun(@times, reshape(reflectance, [r * c, numel(wavelength)]), d65');
    
    XYZ = (CMF' * radiances')';
    XYZ = reshape(XYZ, r,  c, 3);
    XYZ = max(XYZ, 0);
    XYZ = XYZ / max(XYZ(:)); %normalization [0,1]
    
%     if exist('adaptationModel', 'var')
%         targetXYZ = [0.95047; 1.00000; 1.08883]; %d65 white
%         sourceXYZ = [0.0725; 0.0418; 0.1583];
%         M = cbCAT(sourceXYZ, targetXYZ, adaptationModel);
%         
%     else
%         % convert this XYZ representation of the image to the default RGB colour representation sRGB (IEC_61966-2-1),
%         % Forward transformation from 1931 CIE XYZ values to sRGB values (Eqn 6 in
%         % IEC_61966-2-1.pdf).
%         M = [3.2406, -1.5372, -0.4986; ...
%             -0.9689, 1.8758, 0.0414; ...
%             0.0557, -0.2040, 1.0570];
%     end
%     
%     sRGB = (M * XYZ')';
%     
%     % Gamma correction
%     a = 0.055;
%     t = 0.0031308;
%     gamma = 2.4;
%     
%     for i = 1:r * c
%         for j = 1:3
%             if (sRGB(i, j) <= t)
%                 sRGB(i, j) = 12.92 * sRGB(i, j);
%             else
%                 sRGB(i, j) = (1 + a) * sRGB(i, j)^(1 / gamma) - a;
%             end
%         end
%     end
%     
%     % Reshape to recover shape of original input.
%     sRGB = reshape(sRGB, [r, c, 3]);
    
    sRGB = XYZ2sRGB_exgamma(XYZ);

    %clip values to [0,1]
    sRGB = max(sRGB, 0);
    sRGB = min(sRGB, 1);
    
    figure(1);
    imshow(sRGB, 'Border', 'tight');
    title('sRGB image resulting from the MSI image')
    
%     % Plot Chromacity chart
%     x = bsxfun(@(x, s) x./s, XYZ(:, 1), sum(XYZ, 2));
%     y = bsxfun(@(x, s) x./s, XYZ(:, 2), sum(XYZ, 2));
%     figure(2);
%     cieplot();
%     hold on
%     scatter(x, y, 10);
%     hold off
else
    warning('Not implemented')
end

figure(1);
imshow(sRGB);
title('sRGB from MSI data')
% figure(2);
% imshow(sRGB_gamma)
% title('Gamma corrected sRGB from MSI data')


% if ~(isempty(outName))
%     imwrite(sRGB, strcat(outName, '_sRGB', '.tif'), 'tif');
%     
%     %Write to TIFF (in the format that Photoshop can read the L,a,b channels)
%     imwrite(Lab16(:, :, 1), strcat(outName, '_L', '.tif'), 'tif');
%     imwrite(Lab16(:, :, 2)./2+((2^16) / 2), strcat(outName, '_a', '.tif'), 'tif');
%     imwrite(Lab16(:, :, 3)./2+((2^16) / 2), strcat(outName, '_b', '.tif'), 'tif');
% end

%% Makes a 2D histogram style image of wavelength versus intensity of the spectral map

% I = raw2msi(I, 'adjusted');
% 
% figure(2);
% vec = zeros(w, r*c);
% vech = zeros(256, w);
% for i = 1:w
%     vec(i, :) = reshape(squeeze(I(i, :, :)), 1, r*c);
%     h = histogram(vec(i, :).*255, 0:1:255);
%     vech(1:256, i) = h;
% end
% vech = flipud(vech*(2^16 - 1)/max(max(vech)));
% v = figure(3);
% imagesc(vech);
% xlabel('Frequency band')
% ylabel('Image intensity value')
% title('Intensity Histogram')
% colorbar

% outName = fullfile(options.saveOptions.savedir, '9-sRGB', strcat('sRgb_', num2str(id.Group)));
% 
% if ~(isempty(outName))
%     saveas(v, strcat(outName, '_binnedspectrum.jpg'));
% end

end

% ================================================
% *** FUNCTION xyz2srgb
% ***
% *** function [RGB] = xyz2srgb(XYZ)
% *** computes 8-bit sRGB from XYZ
% *** XYZ is n by 3 and in the range 0–1
% *** see also srgb2xyz
function [RGB] = xyz2srgbCustom(XYZ)
if (size(XYZ,2)~=3)
disp('XYZ must be n by 3'); return;
end
M = [3.2406 -1.5372 -0.4986; -0.9689 1.8758 0.0415;
0.0557 -0.2040 1.0570];
RGB = (M*XYZ')';
RGB(RGB>0) = 0;
RGB(RGB<1) = 1;
index = (RGB<=0.0031308);
RGB = RGB + (index).*(12.92*RGB);

RGB = RGB + (1-index).*(1.055*RGB .^ (1/2.4)-0.055);
RGB=ceil(RGB*255);
RGB(RGB<0) = 0;
RGB(RGB>255) = 255;
end
