function dispImage = getDisplayImage(spectralImage, method, channel)
if nargin < 3 
    channel  = 100; 
end 

switch method 
    case 'rgb'
        %[lambda, xFcn, yFcn, zFcn] = colorMatchFcn('CIE_1964');
        [m,n,z] = size(spectralImage);
        lambdaIn = getWavelengths(z, 'raw');    
        [lambdaMatch, xFcn, yFcn, zFcn] = colorMatchFcn('1964_FULL');
        xyz = interp1(lambdaMatch', [xFcn; yFcn; zFcn]', lambdaIn, 'pchip', 0);
        colImage = double(reshape(spectralImage, [m*n,z]));
        %result = arrayfun(@(ROWIDX) naive_bayes(YourArray(ROWIDX,:)), (1:size(YourArray,1)).');
        %colDispImage = applyRowFunc(@(x) squeeze(x) * squeeze(xyz), colImage);
        colDispImage = colImage * squeeze(xyz);
        colDispImage = max(colDispImage, 0);
        colDispImage = colDispImage / max(colDispImage(:));
        dispImage = reshape(colDispImage, [m, n, 3]);
        dispImage = XYZ2sRGB_exgamma(dispImage);
        dispImage = max(dispImage, 0);
        dispImage = min(dispImage, 1);
        dispImage = dispImage.^0.4;
    case 'channel'
        dispImage = rescale(spectralImage(:,:,channel));
    otherwise 
        error('Unsupported method for display image reconstruction');
end 

end 