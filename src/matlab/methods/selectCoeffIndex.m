function [ID] = selectCoeffIndex()
%% SelectCoeffIndex  updates the values of CoeffIndex of ID.mat
% based on the image point which produces minimum average estimation rmse,
% when coefficient values are adjusted to it. 
%
%Usage:
% selectCoeffIndex();

% setSetting('action', 'ReflectanceEstimationPreset');
setSetting('skipLoading', false);
setSetting('showImages', false);
setSetting('saveImages', false);
setSetting('tryReadData', false);
setSetting('pixelValueSelectionMethod', 'adjusted');
setSetting('smoothingMatrixMethod', 'Cor_Sample');
setSetting('noiseType', 'givenSNR');
setSetting('noiseParam', 20);

disp('Searching for reference coefficient...')

%load(generateName('system'), 'wavelength'); % camera system parameters
load(generateName('id'), 'ID'); % image data id and info struct
load(getSetting('matfilein'), 'poiRAWs', 'Spectra');
    
G = findgroups({ID.Sample}, {ID.Type});
for i = 1:max(G)
    rmseMin = 1;
    idxs = find(G == i);
    groupID = ID(idxs);
    [~,coeffIdxs] = unique({groupID.ROI}); %1...N
    coeffIdxs = [groupID(coeffIdxs).Index];
    
    for coeffIndex = coeffIdxs
        avgRmse = 0;
        for k = idxs

            idk = ID(k);
            idk.CoeffIndex = coeffIndex;
            [est, rmse] = estimateReflectance(poiRAWs{k,1}, poiRAWs{k,2}, Spectra(k,:), idk);
            avgRmse = avgRmse + rmse;
        end
        if (avgRmse / length(idxs)) < rmseMin && ~(any(est(:) < 0 ) || any(est(:) > 1 ))
            rmseMin = avgRmse / length(idxs);
            coeffIndexMin = coeffIndex;
        end
    end
    [ID(idxs).CoeffID] = deal(coeffIndexMin);
end

ID = orderfields(ID);
save(fullfile(getSetting('systemdir'), 'ID.mat'), 'ID');
disp('Finished updating CoeffIndex values in the ID file.')

end