load('C:\Users\elena\Google Drive\titech\research\input\saitama_v7_min_region_e\ID.mat');
load('C:\Users\elena\Google Drive\titech\research\input\saitama_v7_min_region_e\data.mat')

dataset = 'saitama_v7_min_region_e';
warning('off')
thresVal = 0.05; 
regionRadius = [];
accTheta = 0.65;
msiN = length(ID);
segmentMaskI = cell(msiN);

basedir = fullfile('..','..','..', 'output', 'general', 'BackgroundRemoval');
saveOptions.saveImages = true;

groups = findgroups([ID.MsiID]);
specimenMasks = cell(max(groups),1);
for g = 1:max(groups)    
    gIdxs = find(groups == g);
    gMembers = ID(gIdxs);
    idd = gMembers(1);
    files = {data([data.MsiID] == idd.MsiID).File};
    [raw, whiteReference, darkReference] = readMSI(files); 
    msi = permute(raw2msi(raw, 'extended'), [2, 3, 1]);
    [m, n, bands] = size(msi);
    bandWeight = 1 / bands; 
    
    %https://www.mathworks.com/help/images/color-based-segmentation-using-k-means-clustering.html
    lab_he = rgb2lab(whiteReference);
%     figure(1);
%     imshow(lab_he);
    ab = lab_he(:,:,2:3);
    ab = im2single(ab);

    nColors = 5;
    % repeat the clustering 3 times to avoid local minima
    pixel_labels = imsegkmeans(ab,nColors,'NumAttempts',3);
%     figure;
%     imshow(pixel_labels,[])
%     title('Image Labeled by Cluster Index');
    maskjj = zeros( m, n);
    bgCounts = zeros(1,nColors);
    for jj = 1:nColors 
%         figure(jj+2);
        maskjj = pixel_labels == jj;
        bgCounts(jj) = sum(sum(maskjj(1:m,1:10))) + sum(sum(maskjj(1:10,1:n)));
%         imshow(maskjj);
    end
    bgCounts 
    bgCounts( bgCounts <= (max(bgCounts)/10)) = 0
    [~,bgChannels] = sort(bgCounts, 'descend');
    bgChannels =  bgChannels(bgCounts(bgChannels) > 0)
    specimenMask =~ismember(pixel_labels, bgChannels);
    
    diam = ceil(min(m,n) * 0.005) 
    se = strel('disk',diam,4);
    specimenMask = imclose(specimenMask,se);
    
    filled = imfill(~specimenMask, 'holes');
    holes = filled & specimenMask;
    bigholes = bwareaopen(holes, 300);
    smallholes = holes & ~bigholes;
    specimenMask =~(~specimenMask | smallholes);

   
    specimenMask = imfill(specimenMask,8, 'holes');
    se = strel('disk',diam*2,4);
    specimenMask = imopen(specimenMask,se);

    specimenMasks{g} = specimenMask;

    cluster1 = whiteReference .* double(specimenMask);
    name = strjoin( {'bgRemoved', num2str(g), idd.Sample, idd.Type}, '_');
    saveOptions.plotName = fullfile(basedir, name);
    plotMontage(whiteReference,cluster1, 2, saveOptions);
    
end
filename = fullfile('..','..','..', 'input', dataset, 'specimenMasks.mat');
save(filename, 'specimenMasks')