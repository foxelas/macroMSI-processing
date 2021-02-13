function resultStruct = prepareMaps2(isFixed, result, mask, id, Iaug, bbox, hasRoi)
   
    roiNames =getSetting('roiNames');
    mapMethods = getSetting('mapMethods');
    
    melMaps = cell(numel(mapMethods), 1);
    hbMaps = cell(numel(mapMethods), 1);
    roiMelMaps = cell(numel(mapMethods), numel(roiNames));
    roiHbMaps = cell(numel(mapMethods), numel(roiNames));
%     mapQualities = zeros(numel(mapMethods), 0);
    
    for i = 1:length(mapMethods)
        close all;
        mmap = result.MelMaps{i,1};
        hbmap = result.HbMaps{i,1};
        [melMaps{i}, hbMaps{i}] = getMap2(mmap, hbmap, mapMethods{i}, mask, id);
        %quality = 0;
        for k = 1:numel(roiNames)
            if hasRoi(k)
                roiMelMap = imcrop(melMaps{i}, bbox{k});
                roiHbMap = imcrop(hbMaps{i}, bbox{k});
                roiMelMaps{i, k} = roiMelMap;
                roiHbMaps{i, k} = roiHbMap;
            end
        end

%         %% GetMap quality
%         %mapQualities(i) = quality / numel(roiNames);
%         if ~any(hasRoi == 0)
%              avgHb = cellfun(@(x) mean(x, 'all'), roiHbMaps(i,:));
%              avgMel = cellfun(@(x) mean(x, 'all'), roiMelMaps(i,:));
% 
%              mapQualities(i) = (avgHb(1) - avgHb(2)) + (avgHb(1) - avgHb(3)) + (avgMel(3) - avgMel(1)) + (avgMel(3) - avgMel(2)); 
%         else
%             mapQualities(i) = nan; 
%         end 
   
    end
    
    resultStruct = struct('IsFixed', isFixed, ... 'MapQualities', {mapQualities}, 
        'MelMaps', {melMaps}, 'HbMaps', {hbMaps}, 'RoiMelMaps', {roiMelMaps}, 'RoiHbMaps', {roiHbMaps}, 'Iaug', Iaug);
end 
