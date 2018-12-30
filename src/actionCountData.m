        %% Count and analyze the contents of the dataset
        [~, idx, ~] = unique(strcat({ID.Csvid}, {ID.T}));
        measuredSpectraCount = length(idx);
        normalFixedCount = 0;
        normalCutCount = 0;
        normalUnfixedCount = 0;
        fixedCount = 0;
        cutCount = 0;
        unfixedCount = 0;
        for i = 1:length(idx)
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
        cancerCount = length(idx) - normalCount;
        cancerUnfixedCount = unfixedCount - normalUnfixedCount;
        fprintf('Breakdown of the dataset:\n')
        fprintf('Total: %d, Unfixed: %d, Fixed: %d, Cut: %d\n', length(idx), unfixedCount, fixedCount, cutCount);
        fprintf('Normal: %d, Cancerous: %d\n\n', normalCount, cancerCount);
        fprintf('Among unfixed data:\n')
        fprintf('Normal: %d, Cancerous: %d\n\n', normalUnfixedCount, unfixedCount-normalUnfixedCount);
        fprintf('Among fixed data:\n')
        fprintf('Normal: %d, Cancerous: %d\n\n', normalFixedCount, fixedCount-normalFixedCount);
        fprintf('Among cut data:\n')
        fprintf('Normal: %d, Cancerous: %d\n\n', normalCutCount, cutCount-normalFixedCount);
        
        %         Breakdown of the dataset:
        %         Unfixed: 43, Fixed: 45, Cut: 39
        %         Normal: 74, Cancerous: 53
        %
        %         Among unfixed data:
        %         Normal: 23, Cancerous: 20
        
        % end of  Count and analyze the contents of the dataset