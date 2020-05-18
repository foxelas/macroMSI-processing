function [ID] = selectCoeffIndex(options)
%% SelectCoeffIndex  updates the values of CoeffIndex of ID.mat
% based on the image point which produces minimum average estimation rmse,
% when coefficient values are adjusted to it. 
%
%Usage:
% selectCoeffIndex('saitama_v2');

% options.action = 'ReflectanceEstimationPreset';
options.skipLoading = false;
options.showImages = false;
options.saveImages = false;
options.tryReadData = false;
options.pixelValueSelectionMethod = 'adjusted';
options.smoothingMatrixMethod = 'Cor_Sample';
options.noiseType = 'givenSNR';
options.noiseParam = 20;

disp('Searching for reference coefficient...')

%load(generateName('system', options), 'wavelength'); % camera system parameters
load(generateName('id', options), 'ID'); % image data id and info struct
load(generateName('matfilein', options), 'poiRAWs', 'Spectra');
    
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
            [est, rmse] = estimateReflectance(poiRAWs{k,1}, poiRAWs{k,2}, Spectra(k,:), idk, options);
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
save(fullfile(options.systemdir, 'ID.mat'), 'ID');
disp('Finished updating CoeffIndex values in the ID file.')

end