function [tables] = exportSimilarityTables(tables, measuredSpectra, adjustedSpectra, alphas, confs)

n = numel(tables);
for i = 1:n
    table1 = tables{i, 1}; 
    [m,z] = size(table1);
    tables{i, 1} = [ table1; [{'Average'}, num2cell(mean(tables{i, 1}{:, 2:end})); ...
                    {'StandardDeviation'}, num2cell(std(tables{i, 1}{:, 2:end}))]];
end

if ~isempty(confs)
    [~,k] = size(confs(1,:));
    for i = 1:n
        tables(i,2:k+1) = deal( confs(i, 1:k));
    end
end 

for j = 1:m
    if tables{1, 1}.AdjGoF(j) < 0.8
        fprintf('**Low result for patch %s\n', tables{1, 1}.Patch{j});
    end
end

savedir = fullfile(getSetting('savedir'), getSetting('experiment'));
save( fullfile(savedir, 'resultTable.mat'), 'tables', 'measuredSpectra', 'adjustedSpectra');
fprintf('Saved .mat values in %s,\n', savedir);

xlsfile = mkNewDir(savedir, 'resultTableFull.xlsx');
for i = 1:n 
    table1 = tables{i,1};
    name = cell(1, size(table1,2));
    name(1:2) = deal(confs(i,:));
    vals =[ name; table1.Properties.VariableNames; table2cell(table1)];
    writecell( vals, xlsfile, 'Sheet', i, 'Range', 'A1','WriteMode', 'overwritesheet');
end 
fprintf('Saved detailed results in %s \n', xlsfile);

xlsfile = mkNewDir(savedir, 'resultTableSummary.xlsx');

name = cell(1, z+1);
name{1} = 'Mean';
vals = [name; 'Configuration', table1.Properties.VariableNames(2:end), 'Multiplier'];
name{1} = 'Mean +/- StD';
vals2 = [name; 'Configuration', table1.Properties.VariableNames(2:end), 'Multiplier'];

for i = 1:n 
    name = {strjoin(confs(i,:), {'_'})};
    meanVals = tables{i, 1}(m+1, 2:end);
    %avg +/ std sheet
    stdVals =  tables{i, 1}(m+2, 2:end);
    meanPlusStds = cellfun(@(x,y) sprintf('%.3f %s %.3f', x, char(177), y) , table2cell(meanVals), table2cell(stdVals), 'un', 0);

    vals = [vals; name, table2cell(meanVals), alphas(i)];
    vals2 = [vals2; name, meanPlusStds, alphas(i)];
end 
writecell( vals, xlsfile, 'Sheet', 1, 'Range', 'A1', 'WriteMode', 'overwritesheet');
writecell( vals2, xlsfile, 'Sheet', 2, 'Range', 'A1', 'WriteMode', 'overwritesheet');
fprintf('Saved summary results in %s \n', xlsfile);



end 
