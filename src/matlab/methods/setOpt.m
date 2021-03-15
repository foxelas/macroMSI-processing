function [] = setOpt(inputSettingsFile)
%     SETOPT configures parameters for running
%
%     Usage:
%     setOpt()
%     setOpt('settings.csv')

if nargin < 1
    inputSettingsFile = '..\..\conf\defaultSettings.csv';
end
tmp = delimread(inputSettingsFile, {', '}, 'raw');
options = struct();
for i = 1:length(tmp.raw)
    
    parameterName = tmp.raw{i, 1};
    rawValue = tmp.raw{i, 2};
    varType = tmp.raw{i, 3};
    
    if isempty(rawValue)
        switch parameterName
            case 'matfilein'
                rawValue = fullfile(options.('systemdir'), 'in.mat');
            case 'matfileinv73'
                rawValue = fullfile(options.('systemdir'), 'in-v73.mat');
            case 'matfileout'
                rawValue = fullfile(options.('savedir'), options.('action'), 'out.mat');
            case 'systemdir'
                rawValue = fullfile('..\..\..\input', options.('dataset'));
            case 'webdir'
                rawValue = fullfile('..\..\..\input', '3_spectra_from_web');
            case 'datadir'
                rawValue = fullfile('..\..\..\..\..\..\..\mspi\', 'saitamav2');
            case 'savedir'
                rawValue = fullfile('..\..\..\output\', options.('dataset'));
        end
    end
    
    if ~isempty(rawValue) && strcmp(parameterName, 'savedir')
        rawValue = fullfile('..\..\..\output\', rawValue);
    end
    
    switch varType
        case 'string'
            value = strrep(rawValue, ' ', ''); %string(rawValue)
        case 'int'
            value = str2num(rawValue);
        case 'double'
            value = str2double(rawValue);
        case 'logical'
            value = strcmp(rawValue, '1');
        case 'doubleArray'
        case 'stringArray'
            value = getStringArray(rawValue);
            
        otherwise
            fprintf('Unsupported type %s for parameter %s.\n', varType, parameterName);
    end
    options.(parameterName) = value;
    eval([parameterName, '=options.', parameterName, ';']);
    
end

options = orderfields(options);

fprintf('Data directory is set to %s.\n', options.datadir);
fprintf('Save directory is set to %s.\n', options.savedir);

clear tmp value parameterName rawValue varType i;
settingsFile = mkNewDir('parameters', 'configuration.mat');
save(settingsFile);
fprintf('Settings loaded from %s and saved in %s.\n', inputSettingsFile, settingsFile);
end

function [arr] = getStringArray(rawValue)
arr = strsplit(rawValue, ' ');
arr = cellfun(@(x) strrep(x, ' ', ''), arr, 'UniformOutput', 0);
end