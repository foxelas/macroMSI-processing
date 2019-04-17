load('D:\temp\Google Drive\titech\research\input\saitama_v7_min_region_e\data.mat')
load('D:\temp\Google Drive\titech\research\input\saitama_v7_min_region_e\ID.mat')
idd = ID(98);
files = {data([data.MsiID] == idd.MsiID).File};    
roiIndexes = find(idd.MsiID == [ID.MsiID]);
rois = length(roiIndexes);
coordinates = [[ID(roiIndexes).Originx]; [ID(roiIndexes).Originy]]';
[segmentMSI, whiteMSI] = readMSI(files);
gg = raw2msi(segmentMSI, 'extended');
[b, m, n] = size(gg);

dataset = 'saitama_v7_min_region_e';
skipLoading = false;
showImages = true;
saveImages = true;
tryReadData = false;

saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
options = struct('tryReadData', tryReadData, 'dataset', dataset, 'action', 'none', ...
    'pixelValueSelectionMethod', 'extended', 'noiseType', 'sameForChannel', 'skipLoading', skipLoading, ...
    'showImages', showImages, 'saveOptions', saveOptions);
options = setOpt(options);

options.saveOptions.BW =true;

%% Show SumLBP
options.saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Texture', strcat(idd.Sample, '_lbp'));

maxScale = 3;
neighbors = 8;
spectralNeighbors = 2;
mapping=getmapping(neighbors,'riu2');

msibands = 9;
scale = 1;
riubins = (neighbors + 2);
sumlbpFeatures = zeros(390,558);

for i = 1:msibands
    lbps = lbp(squeeze(gg(i,:,:)),scale,neighbors,mapping, 'e');
    sumlbpFeatures = sumlbpFeatures + im2double(lbps);
end

lbpFeatures = lbp(rgb2gray(whiteMSI),scale,neighbors,mapping, 'e');
options.saveOptions.BW = false;
plots('lbp', 1, 'Image', lbpFeatures, 'AdditionalImage', sumlbpFeatures, 'SaveOptions', options.saveOptions);
options.saveOptions.BW = true;
plots('lbp', 1, 'Image', lbpFeatures, 'AdditionalImage', sumlbpFeatures, 'SaveOptions', options.saveOptions);

return 