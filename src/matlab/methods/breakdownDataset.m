function [] = breakdownDataset(ID)
%   DATASETBREAKDOWN Counts and analyzes the contents of the dataset
%
%   Input:
%   ID - struct containing dataset information
%
%   Output:
%   void
%
%   Usage:
%   breakdownDataset(ID)

[~, idx, ~] = unique(strcat({ID.SpectrumFile}, {ID.T}));
measuredSpectraCount = length(idx);
normalFixedCount = 0;
normalCutCount = 0;
normalUnfixedCount = 0;
fixedCount = 0;
cutCount = 0;
unfixedCount = 0;
for i = 1:measuredSpectraCount
    if ~(ID(idx(i)).IsFixed)
        unfixedCount = unfixedCount + 1;
        if ID(idx(i)).IsBenign
            normalUnfixedCount = normalUnfixedCount + 1;
        end
    elseif (ID(idx(i)).IsFixed && ~(ID(idx(i)).IsCut))
        fixedCount = fixedCount + 1;
        if ID(idx(i)).IsBenign
            normalFixedCount = normalFixedCount + 1;
        end
    else
        cutCount = cutCount + 1;
        if ID(idx(i)).IsBenign
            normalCutCount = normalCutCount + 1;
        end
    end
end
normalCount = normalUnfixedCount + normalFixedCount + normalCutCount;
cancerCount = measuredSpectraCount - normalCount;
cancerUnfixedCount = unfixedCount - normalUnfixedCount;

fprintf('Breakdown of the dataset:\nTotal: %d, Unfixed: %d, Fixed: %d, Cut: %d\nNormal: %d, Cancerous: %d\n', ...
    length(idx), unfixedCount, fixedCount, cutCount, normalCount, cancerCount);
fprintf('Among unfixed data:\nNormal: %d, Cancerous: %d\n ', normalUnfixedCount, cancerUnfixedCount);
fprintf('Among fixed data:\nNormal: %d, Cancerous: %d\n', normalFixedCount, fixedCount-normalFixedCount);
fprintf('Among cut data:\nNormal: %d, Cancerous: %d\n', normalCutCount, cutCount-normalCutCount);

if ~isempty(getSetting('savedir'))
    fileID = fopen(fullfile(getSetting('savedir'), 'data_count.txt'), 'w');
    
    fprintf(fileID, 'Breakdown of the dataset:\nTotal: %d, Unfixed: %d, Fixed: %d, Cut: %d\nNormal: %d, Cancerous: %d\n', ...
        length(idx), unfixedCount, fixedCount, cutCount, normalCount, cancerCount);
    fprintf(fileID, 'Among unfixed data:\nNormal: %d, Cancerous: %d\n ', normalUnfixedCount, cancerUnfixedCount);
    fprintf(fileID, 'Among fixed data:\nNormal: %d, Cancerous: %d\n', normalFixedCount, cutCount-normalFixedCount);
    fprintf(fileID, 'Among cut data:\nNormal: %d, Cancerous: %d\n', normalCutCount, cutCount-normalCutCount);
    
    fclose(fileID);
end
end
