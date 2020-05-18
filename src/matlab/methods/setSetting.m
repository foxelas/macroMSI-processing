function [] = setSetting(parameter, value)
settingsFile = 'configuration.mat';
m = matfile(settingsFile, 'Writable', true);
if nargin < 2 %write default value
    v = m.options;
    m.(parameter) = v.(parameter);
else
    m.(parameter) = value;
end
end
