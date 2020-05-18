function [value] = getSetting(parameter)
settingsFile = 'configuration.mat';
variableInfo = who('-file', settingsFile);
if ismember(parameter, variableInfo)
    m = matfile(settingsFile);
    value = m.(parameter);
else
    fprintf('Parameter %s does not exist in the configuration file.\n', parameter);
end
end
