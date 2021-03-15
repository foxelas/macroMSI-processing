function [tables] = exportSimilarityTables(tables, measuredSpectra, adjustedSpectra, confs)

n = numel(tables);
for i = 1:n
    tables{i, 1}{25, 1} = {'Average'};
    tables{i, 1}{26, 1} = {'StandardDeviation'};
    for j = 2:7
        tables{i, 1}{25, j} = mean(tables{i, 1}{:, j});
        tables{i, 1}{26, j} = std(tables{i, 1}{:, j});
    end
end

if ~isempty(confs)
    [~,k] = size(confs(1,:));
    for i = 1:n
        tables(i,2:k+1) = deal( confs(i, 1:k));
    end
end 

for j = 1:24
    if tables{1, 1}.AdjGoF(j) < 0.8
        fprintf('**Low result for patch %s\n', tables{1, 1}.Patch{j});
    end
end

savedir = fullfile(getSetting('savedir'), getSetting('experiment'));

save( fullfile(savedir, 'resultTable.mat'), 'tables', 'measuredSpectra', 'adjustedSpectra');
fprintf('Saved .mat values in %s,\n', savedir);

xlsfile = mkNewDir(savedir, 'resultTableFull.xlsx');
for i = 1:n
    writecell( confs(i,:), xlsfile, 'Sheet', i, 'Range', strcat('A1:B', num2str(k)),'WriteMode', 'overwritesheet');
    writetable( tables{i, 1}, xlsfile, 'Sheet', i, 'Range', 'A2');
end

xlsfile = mkNewDir(savedir, 'resultTableSummary.xlsx');
for i = 1:n 
    name = {strjoin(confs(i,:), {'_'})};
    %avg sheet
    writecell( name, xlsfile, 'Sheet', 1, 'Range', strcat('A', num2str(i+2)));
    meanVals = tables{i, 1}(25, 2:end);
    %avg +/ std sheet
    stdVals =  tables{i, 1}(26, 2:end);
    meanPlusStds = cellfun(@(x,y) sprintf('%.3f %s %.3f', x, char(177), y) , table2cell(meanVals), table2cell(stdVals), 'un', 0);
    writecell( name, xlsfile, 'Sheet', 2, 'Range', strcat('A', num2str(i+2)));
    if i == 1 
        writetable( meanVals, xlsfile, 'Sheet', 1, 'Range', strcat('B', num2str(i+1), ':G', num2str(i+2)), 'WriteVariableNames', true);
        writecell( meanVals.Properties.VariableNames, xlsfile, 'Sheet', 2, 'Range', strcat('B', num2str(i+1), ':G', num2str(i+1)));
        writecell( meanPlusStds, xlsfile, 'Sheet', 2, 'Range', strcat('B', num2str(i+2), ':G', num2str(i+2)));
    else 
        writetable( meanVals, xlsfile, 'Sheet', 1, 'Range', strcat('B', num2str(i+2), ':G', num2str(i+2)), 'WriteVariableNames', false);
        writecell( meanPlusStds, xlsfile, 'Sheet', 2, 'Range', strcat('B', num2str(i+2), ':G', num2str(i+2)));
    end 
end 

end 