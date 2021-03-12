function [tables] = exportSimilarityTables(tables, measuredSpectra, adjustedSpectra)

n = numel(tables);
for i = 1:n
    tables{i, 1}{25, 1} = {'Average'};
    tables{i, 1}{26, 1} = {'StandardDeviation'};
    for j = 2:7
        tables{i, 1}{25, j} = mean(tables{i, 1}{:, j});
        tables{i, 1}{26, j} = std(tables{i, 1}{:, j});
    end
end

for j = 1:24
    if tables{1, 1}.AdjGoF(j) < 0.8
        fprintf('**Low result for patch %s\n', tables{1, 1}.Patch{j});
    end
end

savedir = fullfile(getSetting('savedir'), getSetting('saveFolder'), 'resultTable.mat');
save(savedir, 'tables', 'measuredSpectra', 'adjustedSpectra');
fprintf('Saved values in %s,\n', savedir);

for i = 1:n
    writetable(tables{i, 1}, fullfile(getSetting('savedir'), getSetting('saveFolder'), 'resultTable.xlsx'), 'Sheet', i);
end

end 