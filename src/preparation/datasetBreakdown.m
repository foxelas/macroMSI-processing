function [] = datasetBreakdown(ID,options)
%%Count and analyze the contents of the dataset

    [~, idx, ~] = unique(strcat({ID.Csvid}, {ID.T}));
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
            if ID(idx(i)).IsNormal
                normalUnfixedCount = normalUnfixedCount + 1;
            end
        elseif (ID(idx(i)).IsFixed && ~(ID(idx(i)).IsCut))
            fixedCount = fixedCount + 1;
            if ID(idx(i)).IsNormal
                normalFixedCount = normalFixedCount + 1;
            end
        else
            cutCount = cutCount + 1;
            if ID(idx(i)).IsNormal
                normalCutCount = normalCutCount + 1;
            end
        end
    end
    normalCount = normalUnfixedCount + normalFixedCount + normalCutCount;
    cancerCount = measuredSpectraCount - normalCount;
    cancerUnfixedCount = unfixedCount - normalUnfixedCount;

    fprintf('Breakdown of the dataset:\nTotal: %d, Unfixed: %d, Fixed: %d, Cut: %d\nNormal: %d, Cancerous: %d\n',...
            num2str(length(idx)), num2str(unfixedCount),num2str(fixedCount),num2str(cutCount),num2str(normalCount),num2str(cancerCount));
    fprintf('Among unfixed data:\n Normal: %d, Cancerous: %d\n ', num2str(normalUnfixedCount),num2str(cancerUnfixedCount));
    fprintf('Among fixed data:\nNormal: %d, Cancerous: %d\n',  num2str(normalFixedCount), num2str(cutCount-normalFixedCount));
    fprintf('Among cut data:\nNormal: %d, Cancerous: %d\n',  num2str(normalCutCount), num2str(cutCount-normalCutCount));

    if exist('options', 'var') && ~isempty(options.savedir)
        fileID = fopen( fullfile( options.savedir, 'data_count.txt'), 'w');
        
        fprintf(fileID,'Breakdown of the dataset:\nTotal: %d, Unfixed: %d, Fixed: %d, Cut: %d\nNormal: %d, Cancerous: %d\n',...
            num2str(length(idx)), num2str(unfixedCount),num2str(fixedCount),num2str(cutCount),num2str(normalCount),num2str(cancerCount));
        fprintf(fileID,'Among unfixed data:\n Normal: %d, Cancerous: %d\n ', num2str(normalUnfixedCount),num2str(cancerUnfixedCount));
        fprintf(fileID,'Among fixed data:\nNormal: %d, Cancerous: %d\n',  num2str(normalFixedCount), num2str(cutCount-normalFixedCount));
        fprintf(fileID,'Among cut data:\nNormal: %d, Cancerous: %d\n',  num2str(normalCutCount), num2str(cutCount-normalCutCount));
        
        fclose(fileID);
    end
end

