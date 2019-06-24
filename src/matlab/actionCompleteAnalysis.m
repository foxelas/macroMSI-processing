targetIndex = 98;
idd = ID(targetIndex);
files = {data([data.MsiID] == idd.MsiID).File};    
roiIndexes = find(idd.MsiID == [ID.MsiID]);
rois = length(roiIndexes);
coordinates = [[ID(roiIndexes).Originx]; [ID(roiIndexes).Originy]]';

[MSI, RGB] = readMSI(files);
MSI = raw2msi(MSI, 'extended');
[msibands, M, N] = size(MSI);

%% Compute reflectance spectrum 
options.smoothingMatrixMethod = 'Cor_Sample';
options.pixelValueSelectionMethod = 'extended';
options.noiseType = 'spatial olympus';

spect = zeros(rois, length(wavelength));
jj = 0;
for k = roiIndexes
    jj = jj + 1;
    % Give a name
    options.saveOptions.plotName = fullfile(options.saveOptions.savedir, 'CompleteAnalysis', strcat(idd.Sample, '_spect'));
    % Retrieve MSI data
    msi = MSIs{k};
    mask = Masks{k};
    measured = Spectra(k,:);

    [spect(jj,:), ~, ~, ~] = reflectanceEstimation(msi, mask, measured, ID(k), options);                 
    
    if (options.showImages)       
        lines = [measured; spect(jj,:)];    
        lineNames = {'Channel \lambda', 'Measured', 'MSI-Estimated'};
        plots('estimationComparison', 1, lines', idd.Sample, 'Wavelength', wavelength, ...
            'SaveOptions', options.saveOptions, 'LineNames', lineNames);
    end
end  

%% Compute and show texture 
maxScale = 3;
neighbors = 8;
spectralNeighbors = 2;
mapping=getmapping(neighbors,'riu2');
scale = 1;
riubins = (neighbors + 2);

lbpFeats = getLBPFeatures('MMLBP', files, coordinates);
 
sumlbpFeatures = zeros(M-2, N-2);
for i = 1:msibands
    lbps = lbp(squeeze(MSI(i,:,:)),scale,neighbors,mapping, 'e');
    sumlbpFeatures = sumlbpFeatures + im2double(lbps);
end
lbpFeatures = lbp(rgb2gray(RGB),scale,neighbors,mapping, 'e');

if (options.showImages)       
    options.saveOptions.plotName = fullfile(options.saveOptions.savedir, 'CompleteAnalysis', strcat(idd.Sample, '_lbp'));
    plots('LBP', 2, 'Image', lbpFeatures, 'AdditionalImage', sumlbpFeatures, 'SaveOptions', options.saveOptions);
end


%% Plot extracted feature set to python
pause;

%% Visualize Malignancy
currOpt = options;
currOpt.options.saveImages = false;
currOpt.options.showImages = false;
[totalMaskBinary, ~, segmentMaskI] = segmentMSIRegion( files, coordinates, currOpt, 0.75, [], 0.07); %0.8 , 0.08

[M, N, ~] = size(totalMaskBinary);
cancerProb = [0.56506547 0.43493453; 0.39944066 0.60055934; 0.27636779 0.72363221; 0.41575609 0.58424391; 0.57570235 0.42429765; 0.57063083 0.42936917];
isPositive = [1 1 1 1 0 1];
trueLabels = [1 1 1 1 0 0];
options.saveOptions.saveImages = false;
outImgClass = zeros(M,N);

[~, orderedRois] = sort(cellfun(@(x) sum(x(:)), segmentMasks), 'descend');
for i = 1:orderedRois
    [r, c] = find(segmentMaskI{i});
    outImgClass(segmentMaskI{i}) = cancerProb(i,2);
    %imshow(outImg);
    %pause;
end

outImgGroundTruth = zeros(M,N);
for i = 1:rois
    [r, c] = find(segmentMaskI{i});
    outImgGroundTruth(segmentMaskI{i}) = trueLabels(i);
    %imshow(outImg);
end

labels = {isPositive, ~[ID(roiIndexes).IsBenign]};
baseImg = {sRGB, RGB};
overlayImg = {outImgClass, outImgGroundTruth};


if (options.showImages)       
    options.saveOptions.plotName = fullfile(options.saveOptions.savedir, 'CompleteAnalysis', strcat('vis_', idd.Sample));
    plots('visual', 3, labels, '', 'Image', baseImg, 'Overlay', overlayImg, 'Cmap', 'jet', 'Alpha', 0.7, ...
        'Coordinates', coordinates,'saveOptions', options.saveOptions);
end

%%