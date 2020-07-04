% k = 51;
% infile = fullfile(getSetting('systemdir'), 'infiles', strcat('poi_', num2str(k), '.mat'));
% load(infile, 'poiName', 'poiRAW', 'poiSegmentMask', ...
%         'roiSeeds', 'measuredSpectrum', 'poiWhite');
for i = [16, 17]%, 27, 28, 2, 3, 19, 20]
groupId = i;  poiId = groupId2poiId(groupId, ID); %id = ID(78);%fixed data
infile = fullfile(getSetting('systemdir'), 'infiles', strcat('group_', num2str(groupId), '.mat'));
load(infile, 'raw', 'specimenMask');

close all; 
msiType = 'extended'; 
normType = 'divMacbeth';
removeBg = 'true';

setSetting('saveImages', true);
pcComps = getPCMaps(groupId, msiType, normType, removeBg);

if i == 16
    pcComps1 = pcComps; 
end 
if i == 17
    pcComps2 = pcComps; 
end 


end 
%tform1 = getRegistrationTransform(pcComps1, pcComps2, 'surf');
tform2 = getRegistrationTransform(pcComps1, pcComps2, 'regconfig');
recovered = registerImage(pcComps2, pcComps1, tform2);

for i = 1: 3
    pc_fixed = squeeze(pcComps1(i,:,:)); 
    pc_unfixed = squeeze(recovered(i,:,:)); 
    setSetting('plotName', fullfile(getSetting('savedir'), getSetting('pca'), 'ssim', strcat('pc', num2str(i), '.png')));
    plotFunWrapper(1, @plotMontage, pc_unfixed, pc_fixed, 'Registered Unfixed vs Fixed');
    [ssimval,ssimmap] = ssim(pc_unfixed,pc_fixed);
    fprintf('For pc %d: SSIM = %0.5f\n', i, ssimval);
    setSetting('plotName', fullfile(getSetting('savedir'), getSetting('pca'), 'ssim', strcat('pc', num2str(i), '_siim.png')));
    plotFunWrapper(2, @plotSimple, ssimmap, ['Local SSIM Map with Global SSIM Value: ',num2str(ssimval)]);
    pause(0.5);
end 

% [melMap, hbMap] = getMap(raw, 'ding', specimenMask);
% [melMap, hbMap] = getMap(raw, 'vasefi', specimenMask);
%[melMap, hbMap] = getMap(raw, 'ours', specimenMask, id);
%[melMap, hbMap] = getMap(raw, 'diebele', specimenMask, id);

%[melMap, hbMap] = getMap(raw, 'kapsokalyvas', specimenMask, id);

