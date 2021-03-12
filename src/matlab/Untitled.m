
load ref4_scene4.mat;
[m,n,z] = size(reflectances);
slice = reflectances(:,:,17);
figure; imagesc(slice); colormap('gray'); brighten(0.5);

reflectances = reflectances/max(reflectances(:));
 
load illum_6500.mat;
load illum_4000.mat;

lambdaIn = 400:10:720;
[lambdaMatch, xFcn, yFcn, zFcn] = colorMatchFcn('1964_FULL');
xyz = interp1(lambdaMatch', [xFcn; yFcn; zFcn]', lambdaIn, 'pchip', 0);
[solaxSpec, lambdaMatch] = getSolaxSpectra();
illumination = interp1(lambdaMatch, solaxSpec', lambdaIn, 'pchip', 0)';
colImage = reshape(reflectances, m*n, z);
colImage = bsxfun(@times, colImage, illumination');
colXYZ = colImage * squeeze(xyz);
imXYZ = reshape(colXYZ, [m, n, 3]);
imXYZ = max(imXYZ, 0);
imXYZ = imXYZ / max(imXYZ(:));
dispImage = XYZ2sRGB_exgamma(imXYZ);
dispImage = max(dispImage, 0);
dispImage = min(dispImage, 1);
dispImage = dispImage.^0.4;
imshow(dispImage)
   