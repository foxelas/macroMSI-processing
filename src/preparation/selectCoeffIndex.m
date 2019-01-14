function [ID] = selectCoeffIndex(dataset)
%% SelectCoeffIndex  updates the values of CoeffIndex of ID.mat
% based on the image point which produces minimum average estimation rmse,
% when coefficient values are adjusted to it. 
%
%Usage:
% selectCoeffIndex('saitama_v2');

action = lower('ReflectanceEstimationSimple');
skipLoading = false;
showImages = false;
saveImages = false;

saveOptions = struct('saveImages', saveImages, 'saveInHQ', false);
options = struct('tryReadData', false, 'dataset', dataset, 'action', action, ...
'pixelValueSelectionMethod', 'extended', 'noiseType', 'givenSNR', 'snr', 36, 'skipLoading', skipLoading, ...
'showImages', showImages, 'saveOptions', saveOptions);
setup;

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
            g = MSIStruct(k);
            measured = interp1(380:780, measuredSpectrumStruct(k).Spectrum, wavelength, 'nearest');
            idk = ID(k);
            idk.CoeffIndex = coeffIndex;
            [est, rmse] = reflectanceEstimation(g, measured, idk, options);
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
save(fullfile(options.systemdir, 'ID.mat'), 'ID', '-v7.3');
disp('Finished updating CoeffIndex values in the ID file.')

end