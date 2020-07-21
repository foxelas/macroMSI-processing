normType = 'divMacbeth';
removeBg = 'true';
fixedId = 16; 
unfixedId = 17; 

%% Preprocessing for PCA 

% msiType = 'extended'; 
% preprocessingForMaps; 
% 
% %% PCA maps 
% 
% close all; pcCompsFixed = getPCMaps(fixedId, msiType, normType, removeBg, tform2, newDims);
% close all; pcCompsUnfixed = getPCMaps(unfixedId, msiType, normType, removeBg, tform2, newDims);
% 
% %% PCA Similarity 
% close all;
% for i = 1: 3
%     pc_fixed = squeeze(pcCompsFixed(i,:,:)); 
%     pc_unfixed = squeeze(pcCompsUnfixed(i,:,:)); 
%     [ssimval, cc] = getImageSimilarity(pc_unfixed, pc_fixed, strcat('_pc', ...
%         num2str(i)), getSetting('pca'), unfixedId, fixedId, name1, name2);
%     pause(0.5);
% end 

%% Preprocessing for Maps  

% msiType = 'adjusted'; 
% preprocessingForMaps; 

%% Chromophore Absorbance Maps 

id = unfixedId; 
msi = preprocUnfixed;
mask = unfixedMask;
close all; [melMap, hbMap] = getMap(msi, 'ding', mask, id);
close all; [melMap, hbMap] = getMap(msi, 'vasefi', mask, id);
close all; [melMap, hbMap] = getMap(msi, 'ours', mask, ID(id));

id = fixedId; 
msi = preprocFixed;
mask = fixedMask;
close all; [melMap, hbMap] = getMap(msi, 'ding', mask, id);
close all; [melMap, hbMap] = getMap(msi, 'vasefi', mask, id);
close all; [melMap, hbMap] = getMap(msi, 'ours', mask, ID(id));


%% Chromophore Index Maps 

% id = unfixedId; 
% msi = preprocUnfixed;
% mask = unfixedMask;
% close all; [melMap, hbMap] = getMap(msi, 'diebele', mask, id);
% close all; [melMap, hbMap] = getMap(msi, 'kapsokalyvas', mask, id);
% close all; [melMap, hbMap] = getMap(msi, 'kuzmina', mask, id);
% 
% id = fixedId; 
% msi = preprocFixed;
% mask = fixedMask;
% close all; [melMap, hbMap] = getMap(msi, 'diebele', mask, id);
% close all; [melMap, hbMap] = getMap(msi, 'kapsokalyvas', mask, id);
% close all; [melMap, hbMap] = getMap(msi, 'kuzmina', mask, id);





