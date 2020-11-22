function nameList = combineNameLists(namesRows, namesColumns)
nameList = cell(numel(namesRows)*numel(namesColumns), 1);
for i = 1:numel(namesRows)
    for j = 1:numel(namesColumns)
        nameList{(i - 1)*numel(namesColumns)+j} = strcat(namesColumns{j}, '_', namesRows{i});
    end
end

end 