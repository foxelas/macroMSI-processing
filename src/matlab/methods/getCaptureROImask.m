function mask = getCaptureROImask(m,n)
%%GETCAPTUREROIMASK returns the mask for the ROI that was selected during
%%catpure.
%   Usage:
%   [m,n,w] = spectralData; 
%   mask = getCaptureROIMask(m,n);

mask = zeros(1376, 1024);
maxSize = [1376, 1024];
if isequal([m,n], maxSize)
    mask = ones(maxSize);
elseif belongsInAvailableROIs(m,n)
    v = floor((maxSize - [m,n]) / 2);  
    mask(v(1):m, v(2):n) = 1;    
%     mask = mask(floor(m/2)-floor(x/2):floor(m/2)+floor(x/2)-1, floor(n/2)-floor(y/2):floor(n/2)+floor(y/2)-1);
elseif isequal([m,n], [1088, 982])
    mask = ones([1088, 982]);
else 
    error('Unsupported capture ROI.');
end 


end

function roiBelongs = belongsInAvailableROIs(m,n)
    roiBelongs = (isequal([m, n], [688, 512]) || isequal([m,n], [458, 341]));
end 