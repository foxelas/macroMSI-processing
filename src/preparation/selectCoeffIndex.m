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
options.pixelValueSelectionMethod = 'extended';
options.smoothingMatrixMethod = 'Cor_Sample';

% out = matfile(options.outName, 'Writable', true);
% readData;
load(fullfile(options.systemdir, 'system.mat'), 'wavelength'); % camera system parameters
% load(fullfile(options.systemdir, 'data.mat')); % image data
load( fullfile(options.systemdir, 'ID.mat'), 'ID'); % image data id and info struct
in = matfile(generateName(options, 'matfilein'));
    
MSIStruct = in.MSIStruct;
measuredSpectrumStruct = in.MeasuredSpectrumStruct;
    
    
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
            measured = measuredSpectrumStruct(k).Spectrum; 
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