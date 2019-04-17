load('D:\temp\Google Drive\titech\research\input\saitama_v7_min_region_e\data.mat')
load('D:\temp\Google Drive\titech\research\input\saitama_v7_min_region_e\ID.mat')
idd = ID(85);
files = {data([data.MsiID] == idd.MsiID).File};    
roiIndexes = find(idd.MsiID == [ID.MsiID]);
rois = length(roiIndexes);
coordinates = [[ID(roiIndexes).Originx]; [ID(roiIndexes).Originy]]';
[segmentMSI, viewImg] = readMSI(files);

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

options.saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Texture', strcat(idd.Sample, '_lbp'));
        
currOpt = options;
currOpt.options.saveImages = false;
[totalMaskBinary, totalMaskColor, segmentMaskI] = segmentedRegions( files, coordinates(:,:), currOpt, 0.75, [], 0.07, []); %0.8 , 0.08

[M, N, ~] = size(totalMaskBinary);
cancerProb = [0.56506547 0.43493453; 0.39944066 0.60055934; 0.27636779 0.72363221; 0.41575609 0.58424391; 0.57570235 0.42429765; 0.57063083 0.42936917];
isPositive = [1 1 1 1 0 1];
trueLabels = [1 1 1 1 0 0];
options.saveOptions.saveImages = false;
outImg = zeros(M,N);
for i = 1:rois
    [r, c] = find(segmentMaskI{i});
    outImg(segmentMaskI{i}) = cancerProb(i,2);
    %imshow(outImg);
    %pause;
end
options.saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Visualization', strcat('pred_', idd.Sample));
plots('visual', 2, isPositive, 'Malignancy Probability', 'Image', viewImg, 'Overlay', outImg, 'Cmap', 'jet', 'Alpha', 0.7, ...
    'Coordinates', coordinates,'saveOptions', options.saveOptions);

outImg = zeros(M,N);
for i = 1:rois
    [r, c] = find(segmentMaskI{i});
    outImg(segmentMaskI{i}) = trueLabels(i);
    %imshow(outImg);
end
options.saveOptions.plotName = fullfile(options.saveOptions.savedir, 'Visualization', strcat('gt_', idd.Sample));
plots('visual', 3, ~[ID(roiIndexes).IsBenign], 'Ground Truth', 'Image', viewImg, 'Overlay', outImg, 'Cmap', 'jet', 'Alpha', 0.7, ...
     'Coordinates', coordinates, 'saveOptions', options.saveOptions);


% options.pixelValueSelectionMethod = 'extended';
% options.smoothingMatrixMethod = 'Cor_Sample';
% options.noiseType = 'sameForChannel 0.0001';  
% estimated = reflectanceEstimation(segmentMSI, [], [], idd, options);
% [W, M, N] = size(estimated);
% estimated_2D = reshape(estimated, W, M*N)';
% save(strcat(options.saveOptions.savedir, '\img_spectra.mat'),'estimated_2D');
% load(strcat(options.saveOptions.savedir, '\predictions.mat'))
% predictions = reshape(predictions, M, N);
% figure(1);imshow(predictions)
% scores = reshape(scores(:,1), M, N);
% figure(2);imshow(scores)



