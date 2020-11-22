function resultStruct = prepareMaps(isFixed, msi, mask, id, Iaug, bbox, hasRoi)
   
    roiNames =getSetting('roiNames');
    mapMethods = getSetting('mapMethods');
    
    melMaps = cell(numel(mapMethods), 1);
    hbMaps = cell(numel(mapMethods), 1);
    roiMelMaps = cell(numel(mapMethods), numel(roiNames));
    roiHbMaps = cell(numel(mapMethods), numel(roiNames));
    mapQualities = zeros(numel(mapMethods), 0);
    
    for i = 1:length(mapMethods)
        close all;
        [melMaps{i}, hbMaps{i}] = getMap(msi, mapMethods{i}, mask, id);
        %quality = 0;
        for k = 1:numel(roiNames)
            if hasRoi(k)
                roiMelMap = imcrop(melMaps{i}, bbox{k});
                roiHbMap = imcrop(hbMaps{i}, bbox{k});
                roiMelMaps{i, k} = roiMelMap;
                roiHbMaps{i, k} = roiHbMap;

    %             if strcmp(roiNames{i}, 'Mel')
    %                 meanMel = mean(roiMelMap(:)); 
    %                 %quality = quality + sum(roiMelMap(:) > (1 + 0.2)*mean(roiMelMap(:))) / length(roiMelMap(:));
    %             elseif strcmp(roiNames{i}, 'Hb')
    %                 meanHb = [ mean(roiHbMap(:))];
    %                 %quality = quality + sum(roiHbMap(:) > (1 + 0.1)*mean(roiHbMap(:))) / length(roiMelMap(:));
    %             else
    %                 meanNorm = mean(roiMelMap(:)); 
    %                 %quality = quality + 1 - sum(roiMelMap(:) > (1 + 0.2)*mean(roiMelMap(:))) / length(roiMelMap(:)) ...
    %                 %    + 1 - sum(roiHbMap(:) > (1 + 0.1)*mean(roiHbMap(:))) / length(roiMelMap(:));
    %             end
            end
        end

        %% GetMap quality
        %mapQualities(i) = quality / numel(roiNames);
        if ~any(hasRoi == 0)
             avgHb = cellfun(@(x) mean(x, 'all'), roiHbMaps(i,:));
             avgMel = cellfun(@(x) mean(x, 'all'), roiMelMaps(i,:));

             mapQualities(i) = (avgHb(1) - avgHb(2)) + (avgHb(1) - avgHb(3)) + (avgMel(3) - avgMel(1)) + (avgMel(3) - avgMel(2)); 
        else
            mapQualities(i) = nan; 
        end 
   
    end
    
    resultStruct = struct('IsFixed', isFixed, 'MapQualities', {mapQualities}, 'MelMaps', {melMaps}, ...
        'HbMaps', {hbMaps}, 'RoiMelMaps', {roiMelMaps}, 'RoiHbMaps', {roiHbMaps}, 'Iaug', Iaug);
end 
