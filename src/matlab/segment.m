load('C:\Users\elena\Google Drive\titech\research\input\saitama_v7_min_region_e\ID.mat');
load('C:\Users\elena\Google Drive\titech\research\input\saitama_v7_min_region_e\data.mat')
load('C:\Users\elena\Google Drive\titech\research\input\saitama_v7_min_region_e\specimenMasks.mat')

warning('off')
thresVal = 0.05; 
regionRadius = [];
accTheta = 0.65;
msiN = length(ID);
segmentMaskI = cell(msiN);

dataset = 'saitama_v7_min_region_e';
basedir = fullfile('..','..','..', 'output', dataset);
saveOptions.saveImages = true;

%add grou[s to ID 
groups = findgroups([ID.MsiID]);
for k = 6:msiN   
    idd = ID(k);
    files = {data([data.MsiID] == idd.MsiID).File};
    [raw, whiteReference, darkReference] = readMSI(files); 
    msi = permute(raw2msi(raw, 'extended'), [2, 3, 1]);
    [m, n, bands] = size(msi);
    bandWeight = 1 / bands; 
        
    maskAgreement = zeros(m, n); % Final mask for the region 
    x = idd.Originx;
    y = idd.Originy;
    for i = 1:bands  
        I2 = squeeze(msi(:,:,i)) .* specimenMasks{groups(k)};
        [~, maskTmp]= regionGrowing(I2, [y, x], thresVal, regionRadius);
%         figure(i);
%         montage({I2, double(maskTmp)});
        maskAgreement = maskAgreement + bandWeight * maskTmp; 
    end

    mask = maskAgreement > 0.2; %accTheta; 
    mask = imfill(mask, 'holes');
    diam = ceil(min(m,n) * 0.005) 
    se = strel('disk',diam,4);
    mask = imclose(mask, se);
    if sum(mask(:)) < 9 %if the region is too small, add the neigborhood of the region seed pixel
        patchX = (x-2):(x+2);
        patchY = (y-2):(y+2);
        mask(patchY, patchX) = 1;
    end
    segmentMaskI{k} = mask;
    saveOptions.saveImages = true;
    saveOptions.saveInHQ = false;
    if (idd.IsBenign)
        mal = 'normal';
    else 
        mal = 'cancer';
    end
    saveOptions.plotName = strcat('..\..\..\output\segments\', 'segment_', num2str(k), '_', idd.Sample, '_', mal, '.png');
    plots('segmentation', 10, 'Image', whiteReference, 'Overlay', maskAgreement, ...
            'AdditionalImage',  mask, 'Coordinates', [x,y], 'SaveOptions', saveOptions);
  
end
%filename = fullfile('..','..','..', 'input', 'saitama_v7_min_region_e', 'specimenMasks.mat');
%save(filename, 'specimenMasks')