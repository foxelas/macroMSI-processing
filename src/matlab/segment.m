function segmentMask = segment(raw, idd, whiteReference, specimenMask, options, regionRadius, thresVal, accTheta)
%load('C:\Users\elena\Google Drive\titech\research\input\saitama_v7_min_region_e\ID.mat');
%load('C:\Users\elena\Google Drive\titech\research\input\saitama_v7_min_region_e\data.mat')
%load('C:\Users\elena\Google Drive\titech\research\input\saitama_v7_min_region_e\specimenMasks.mat')
%dataset = 'saitama_v8_min_region_bright';

if (nargin < 6)
    regionRadius = [];
end
if (nargin < 7)
    thresVal = 0.05; 
end
if (nargin < 8)
    accTheta = 1/3;
end 

isSmallMask = ~isempty(thresVal);
 
saveOptions = options.saveOptions;
saveOptions.saveImages = true;
saveOptions.saveInHQ = false;

msi = permute(raw2msi(raw, 'extended'), [2, 3, 1]);
[m, n, bands] = size(msi);
bandWeight = 1 / bands; 

maskAgreement = zeros(m, n); % Final mask for the region 
x = idd.Originx;
y = idd.Originy;
diam = ceil(min(m,n) * 0.01);
se = strel('disk',diam,4);

for i = 1:bands  
    I2 = squeeze(msi(:,:,i)) .* specimenMask;
    [~, maskTmp]= regionGrowing(I2, [y, x], thresVal, regionRadius);
    maskTmp = imdilate(maskTmp,se);
    maskTmp = imclose(maskTmp, se);

%         figure(i); %
%         montage({I2, double(maskTmp)}); %
    maskAgreement = maskAgreement + bandWeight * maskTmp; 
end

mask = maskAgreement > accTheta; 
mask = ~bwareaopen(~mask, ceil(m*n/100), 8);
mask = imclose(mask, se);
mask = imopen(mask, se);

if sum(mask(:)) < 9 %if the region is too small, add the neigborhood of the region seed pixel
    patchX = (x-5):(x+5);
    patchY = (y-5):(y+5);
    mask(patchY, patchX) = 1;
end

segmentMask = mask .* specimenMask;

if (idd.IsBenign)
    mal = 'benign';
else 
    mal = 'malignant';
end

if (isSmallMask)
    folder =  outputFolderMap('segmentsForFeatureExtraction');
else
    folder = outputFolderMap('segments');
end
saveOptions.plotName = fullfile(saveOptions.savedir, folder, strcat( 'segment_', num2str(idd.Index), '_', idd.Sample, '_', mal));
baseImage = whiteReference .* specimenMask;
plotSegmentation(baseImage, maskAgreement, baseImage .* mask, [x,y], 1, saveOptions);

end