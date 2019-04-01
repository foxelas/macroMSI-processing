load('C:\Users\elena\Google Drive\titech\research\input\saitama_v5_min_region\data.mat')
load('C:\Users\elena\Google Drive\titech\research\input\saitama_v5_min_region\ID.mat')
load('temp.mat')
idd = ID(369);
files = {data([data.MsiID] == idd.MsiID).File};    
roiIndexes = find(idd.MsiID == [ID.MsiID]);
rois = length(roiIndexes);
coordinates = [[ID(roiIndexes).Originx]; [ID(roiIndexes).Originy]]';
    
[segmentMSI, whiteMSI] = readMSI(files);
viewImg = whiteMSI;
figure(3); 
imshow(viewImg)

colors = jet(size(coordinates,1));
hold on
for i = 1:size(coordinates, 1)
    x = coordinates(i,1);
    y = coordinates(i,2);
    plot(x, y, '*', 'LineWidth', 2, 'MarkerSize', 5, 'MarkerFaceColor', colors(i,:));
end
hold off
        
dataset = 'saitama_v5_min_region';
skipLoading = false;
showImages = true;
saveImages = false;
tryReadData = false;

saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
options = struct('tryReadData', tryReadData, 'dataset', dataset, 'action', 'none', ...
    'pixelValueSelectionMethod', 'extended', 'noiseType', 'sameForChannel', 'skipLoading', skipLoading, ...
    'showImages', showImages, 'saveOptions', saveOptions);
options = setOpt(options);

%[totalMaskBinary, totalMaskColor, segmentMaskI] = segmentedRegions( files, coordinates, options, 0.6, [], [], 0.08, []);
[M, N, ~] = size(viewImg);
cancerProb = [ 0.6711745  0.3288255 ; 
 0.3576419  0.6423581;
 0.27323704 0.72676296;
 0.23336374 0.76663626;
 0.91281083 0.08718917;
 0.26809716 0.73190284];

isPositive = [0 1 1 1 0 1];
trueLabels = [1 1 1 1 0 0];
 
% figure(4);
% hold on
% outImg = zeros(M,N,3);
% for i = 1:rois
%     outImg = outImg +  cat(3, double(segmentMaskI{i}), double(segmentMaskI{i}),double(segmentMaskI{i}));
%     imshow(outImg);
%     pause(0.5)
% end
% hold off

figure(5);
outImg = zeros(M,N);
for i = 1:rois
    [r, c] = find(segmentMaskI{i});
    outImg(segmentMaskI{i}) = cancerProb(i,2);
    imshow(outImg);
end
plots('visual', 1, 'Image', viewImg, 'Overlay', outImg, 'Cmap', 'jet', 'Alpha', 0.7, 'Title', 'Malignancy Probability');

outImg = zeros(M,N);
for i = 1:rois
    [r, c] = find(segmentMaskI{i});
    outImg(segmentMaskI{i}) = isPositive(i);
    imshow(outImg);
end
plots('visual', 2, 'Image', viewImg, 'Overlay', outImg, 'Cmap', 'jet', 'Alpha', 0.7,  'Title', 'Classification');

outImg = zeros(M,N);
for i = 1:rois
    [r, c] = find(segmentMaskI{i});
    outImg(segmentMaskI{i}) = trueLabels(i);
    imshow(outImg);
end
plots('visual', 3, 'Image', viewImg, 'Overlay', outImg, 'Cmap', 'jet', 'Alpha', 0.7,  'Title', 'Ground Truth');


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

% [M, N, ~] = size(whiteMSI);
% estimated = zeros(M,N,81);
% step = 5;
% for i = step:(M-step)
%     for j = step:(N-step)
%         estimated(i,j,:) = reflectanceEstimation(squeeze(segmentMSI(:, i-step:i+step, j-step:j+step, :)), [], [], idd, options);
% 
%     end
% end

% currentOptions.showImages = false; 
% 
% roiIndexes = find(idd.RgbID == [ID.RgbID]);
% rois = length(roiIndexes);
% roiNames = cell(rois,1);
% roiPositions = zeros(rois,2);
% roiMalignancy = zeros(rois,1);
% for i=1:rois
%     roiNames{i} = strcat(ID(roiIndexes(i)).SpectrumFile,'_' ,ID(roiIndexes(i)).T);
%     roiPositions(i, :) = [ID(roiIndexes(i)).Originx, ID(roiIndexes(i)).Originy];
%     roiMalignancy(i) = ~ID(roiIndexes(i)).IsBenign;
%     
%     [segmentMSI, segmentWhite, segmentDark, segmentMask, segmentMaskI] = segmentMSIRegion(files, roiPositions(i,:), currentOptions);
%     plots('segmentation', 2, 'Image', viewImg + segmentMaskI{1}, 'Coordinates', roiPositions(i,:));
%     figure(3);
%     imoverlay(rgb2gray(viewImg), double(segmentMaskI{1}), [], [], 'parula')
%     
%     pause; 
% 
% end

