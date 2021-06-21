function [] = setSetting(parameter, value)
%     SETSETTING sets a parameter according to a value or by default
%
%     Usage:
%     setSetting('savedir', 'out\out')
%     setSetting('savedir')

settingsFile = mkNewDir('parameters', 'configuration.mat');
m = matfile(settingsFile, 'Writable', true);
if nargin < 2 %write default value
    v = m.options;
    m.(parameter) = v.(parameter);
    value = v.(parameter);
else
    m.(parameter) = value;
end
notifySetting(parameter, value);
end
