function specimenMask = removeBackground(whiteReference, idd, colorLevelsForKMeans, attemptsForKMeans, bigHoleCoefficient, closingCoefficient, openingCoefficient)

%load( fullfile('..','..','..', 'input\saitama_v7_min_region_e\ID.mat'));
%load( fullfile('..','..','..', 'input\saitama_v7_min_region_e\data.mat'));
%dataset = 'saitama_v7_min_region_e';

if (nargin < 4)
    colorLevelsForKMeans = 6;
end
if (nargin < 5)
    attemptsForKMeans = 3;
end
if (nargin < 6)
    layerSelectionThreshold = 0.1;
end
if (nargin < 7)
    bigHoleCoefficient = 100;
end
if (nargin < 8)
    closingCoefficient = 2;
end
if (nargin < 9)
    openingCoefficient = 5;
end

[m, n, ~] = size(whiteReference);

%https://www.mathworks.com/help/images/color-based-segmentation-using-k-means-clustering.html
lab_he = rgb2lab(whiteReference);
%     figure(1);
%     imshow(lab_he);
ab = lab_he(:, :, 2:3);
ab = im2single(ab);

% repeat the clustering 3 times to avoid local minima
pixel_labels = imsegkmeans(ab, colorLevelsForKMeans, 'NumAttempts', attemptsForKMeans);
%     figure(13);
%     imshow(pixel_labels,[])
%     title('Image Labeled by Cluster Index');

bgCounts = zeros(1, colorLevelsForKMeans);
for jj = 1:colorLevelsForKMeans
    %         figure(jj+2); %%
    maskjj = pixel_labels == jj;
    bgCounts(jj) = sum(sum(maskjj(1:m, 1:10))) + sum(sum(maskjj(1:10, 1:n)));
    %         imshow(maskjj); %%
end
%bgCounts
bgCounts(bgCounts <= (max(bgCounts) * layerSelectionThreshold)) = 0; %%
[~, bgChannels] = sort(bgCounts, 'descend');
bgChannels = bgChannels(bgCounts(bgChannels) > 0); %%
specimenMask = ~ismember(pixel_labels, bgChannels);

diam = ceil(min(m, n)*0.005); %%
specimenMask = imclose(specimenMask, strel('disk', diam*closingCoefficient, 4));

specimenMask = closeSmallHolesInTheBackground(specimenMask, diam*bigHoleCoefficient);

specimenMask = imfill(specimenMask, 8, 'holes');
specimenMask = imopen(specimenMask, strel('disk', diam*openingCoefficient, 4));
specimenMask = bwareaopen(specimenMask, ceil(m*n/500), 8);
%specimenMask = imopen(specimenMask,strel('disk',diam * openingCoefficient,4));

cluster1 = whiteReference .* double(specimenMask);
name = strjoin({'bgRemoved', num2str(idd.Group), idd.Sample, idd.Type}, '_');
 setSetting( 'plotName', fullfile(getSetting('savedir'), getSetting('backgroundRemoval'), name));
plotMontage(whiteReference, cluster1, 1);

end

function imageWithoutSmallHoles = closeSmallHolesInTheBackground(image, bigHoleDiameter)
filled = imfill(~image, 'holes');
holes = filled & image;
bigholes = bwareaopen(holes, bigHoleDiameter);
smallholes = holes & ~bigholes;
imageWithoutSmallHoles = ~(~image | smallholes);
end
