function [refXYZ] = referenceXYZ(pixelValueSelectionMethod, illumination)

wavelength = [450, 465, 505, 525, 575, 605, 630];
flist = dir('D:\temp\mspi\saitamav2\white\unfixed\left');
g = readMSI({flist(12:18).name}, 826, 2445, 30, 30, wavelength);
figure;
imshow(squeeze(g(1, :, :, :)));
gg = valueSelect(g, pixelValueSelectionMethod);
refWhite = mean(reshape(gg, [30 * 30, 7]), 1);

% load color matching functions
[lambdaMatch, xFcn, yFcn, zFcn] = colorMatchFcn('1931_FULL');

% Interpolate the input wavelength in the color matching functions. 7x3
xyzbar = interp1(lambdaMatch', [xFcn; yFcn; zFcn]', wavelength', 'pchip', 0);

% Illumination
if size(illumination, 1) == 401
    lambda = 380:780;
elseif size(illumination, 1) == 81
    lambda = 380:5:780;
else
    lambda = 300:800;
end

illum = interp1(lambda, illumination, wavelength', 'nearest');

radiance = refWhite * illum;

refXYZ = (xyzbar' * radiance');
%refXYZ = refXYZ/max(refXYZ); %normalization [0,1]
end
