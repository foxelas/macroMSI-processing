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

load(fullfile(options.systemdir, 'system.mat'), 'wavelength'); % camera system parameters
load(fullfile(options.systemdir, 'ID.mat'), 'ID'); % image data id and info struct
load(generateName(options, 'matfilein'), 'MSIs', 'Masks', 'Spectra');
    
G = findgroups([ID.Sample], [ID.Type]);
for i = 1:max(G)
    rmseMin = 1;
    idxs = find(G == i);
    groupID = ID(idxs);
    [~,coeffIdxs] = unique({groupID.IMG}); %1...N
    coeffIdxs = [groupID(coeffIdxs).Index];
    
    for coeffIndex = coeffIdxs
        avgRmse = 0;
        for k = idxs

            idk = ID(k);
            idk.CoeffIndex = coeffIndex;
            [est, rmse] = reflectanceEstimation(MSIs{k}, Masks{k}, Spectra(k,:), idk, options);
            avgRmse = avgRmse + rmse;
        end
        if (avgRmse / length(idxs)) < rmseMin && ~(any(est(:) < 0 ) || any(est(:) > 1 ))
            rmseMin = avgRmse / length(idxs);
            coeffIndexMin = coeffIndex;
        end
    end
    [ID(idxs).CoeffIndex] = deal(coeffIndexMin);
end

ID = orderfields(ID);
save(fullfile(options.systemdir, 'ID.mat'), 'ID');
disp('Finished updating CoeffIndex values in the ID file.')

end