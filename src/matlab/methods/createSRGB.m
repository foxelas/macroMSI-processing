function [sRGB, Lab16] = createSRGB(I, method, id, adaptationModel, mask)
%   CREATESGB generates the sRGB based on the MSI reflectances
%
%   Usage: sRGB = createSRGB( I, 'medium', id, 'vonKries', mask);


wavelength = 380:5:780;
[~, r, c, ~] = size(I);
if nargin < 6
    mask = ones(r, c);
end
% load color matching functions
[lambdaMatch, xFcn, yFcn, zFcn] = colorMatchFcn('1931_FULL');
[~, wavelengthIdx, ~] = intersect(lambdaMatch, wavelength);
CMF = [xFcn(wavelengthIdx); yFcn(wavelengthIdx); zFcn(wavelengthIdx)]';

% Illumination
[lambdaMatch, d65full] = illuminant('d65');
[~, wavelengthIdx, ~] = intersect(lambdaMatch, wavelength);
d65 = d65full(wavelengthIdx);

setSetting('smoothingMatrixMethod', 'Cor_All');
setSetting('pixelValueSelectionMethod', 'extended');
setSetting('noiseType', 'spatiospectralolympus');
setSetting('rho', 0.6);
setSetting('windowDim', 3);
setSetting('noiseParam', 0.000001); % 0.0001;

% setSetting('noiseType', 'fromOlympus');
% setSetting('noiseParam', 0.001);

reflectance = estimateReflectance(I, [], [], id);
reflectance = permute(reflectance, [2, 3, 1]);

XYZideal = bsxfun(@times, CMF, d65);
%Compute the normalisation factor k
reflectance = reshape(reflectance, r*c, numel(wavelength));
XYZobj = [reflectance * XYZideal(:, 1), reflectance * XYZideal(:, 2), reflectance * XYZideal(:, 3)];
k = 1 / sum(XYZideal(:, 2));
XYZobj_col = k * XYZobj;

if strcmp(method, 'original')
    
    sRGB_col = xyz2srgbCustom(XYZobj_col);
    sRGB = reshape(sRGB_col, r, c, 3);
    
elseif strcmp(method, 'medium')
    sourceXYZ = [0.2; 0.2; 0.15]; %[0.1876; 0.1928; 0.1756];
    targetXYZ = [1.09846607; 1.00000000; 0.35582280]; %d65 white
    Madapt = cbCAT(sourceXYZ, targetXYZ, adaptationModel);
    
    sRGBadapt = (Madapt * XYZobj_col')';
    sRGB = reshape(sRGBadapt, r, c, 3);
    
else
    warning('Not implemented')
end

figure(1);
sRGB = sRGB .* mask;
imshow(sRGB);
title('sRGB from MSI data')

%{
%     % Plot Chromacity chart
%     x = bsxfun(@(x, s) x./s, XYZ(:, 1), sum(XYZ, 2));
%     y = bsxfun(@(x, s) x./s, XYZ(:, 2), sum(XYZ, 2));
%     figure(2);
%     cieplot();
%     hold on
%     scatter(x, y, 10);
%     hold off


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

% outName = fullfile(getSetting('savedir'), getSetting('sRGB'), strcat('sRgb_', num2str(id.Group)));
%
% if ~(isempty(outName))
%     saveas(v, strcat(outName, '_binnedspectrum.jpg'));
% end
%}
end


function xy = XYZ2xy(xyz)
%xy = XYZ2xy(xyz)
% Converts CIE XYZ to xy chromaticity.

X = xyz(1);
Y = xyz(2);
s = sum(xyz);
xy = [X / s; Y / s];

end

function xyz = xy2XYZ(xy, Y)
%xyz = xy2XYZ(xy,Y)
% Converts xyY chromaticity to CIE XYZ.
x = xy(1);
y = xy(2);
xyz = [Y / y * x; Y; Y / y * (1 - x - y)];
end