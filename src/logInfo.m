function [outputLog] = logInfo(options, runTime)
%LOGINFO produces a log of the current action run

%% Make a struct without nested structs
saveOptions = options.saveOptions;
f = fieldnames(saveOptions);
for i = 1:length(f)
    options.(f{i}) = saveOptions.(f{i});
end
options = rmfield(options, {'saveOptions'});

%% Produce output log
optionsTab = struct2table(options);
outputLog = '-------------------------------------------\n';
outputLog = strjoin([outputLog, 'Run on ', string(datetime), '\n'], ' ');
outputLog = strjoin([outputLog, string(runTime), 's elapsed.', '\n'], ' ');
for i = 1:width(optionsTab)
    outputLog = strjoin([outputLog, optionsTab.Properties.VariableNames(i), strrep(strrep(string(optionsTab(1, i).Variables), '..', ''), '\', ' '), '\n'], '     ');
end
outputLog = strjoin([outputLog, '-------------------------------------------\n\n\n\n\n\n\n']);
fileID = fopen('..\logs\log.txt', 'a');
fprintf(fileID, outputLog);
fclose(fileID);

end
