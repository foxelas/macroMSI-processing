function [] = setOpt(inputSettingsFile)

if nargin < 1 
    inputSettingsFile = '..\..\conf\defaultSettings.csv';
end
tmp = delimread(inputSettingsFile, ', ', 'raw');
options = struct();
for i = 1:length(tmp.raw)
    
    parameterName = tmp.raw{i,1};
    rawValue = tmp.raw{i,2};
    varType = tmp.raw{i,3};
    
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
            case 'datadir'
                rawValue = fullfile('.\..\..\..\..\..\mspi\', options.('dataset'));
            case 'savedir'
                rawValue = fullfile('..\..\..\output\', options.('dataset'));       
        end
    end 
    
    switch varType
        case 'string'
            value = rawValue; %string(rawValue)
        case 'int'
            value = str2num(rawValue);
        case 'double'
            value = str2double(rawValue);
        case 'logical'
        	value = strcmp(rawValue, '1');
        otherwise 
            fprintf('Unsupported type %s for parameter %s.\n', varType, parameterName);
    end 
    options.(parameterName) = value;
    eval([parameterName '=options.' parameterName ';']);

end 

options = orderfields(options);

fprintf('Data directory is set to %s.\n', options.datadir);
fprintf('Save directory is set to %s.\n', options.savedir);

clear tmp value parameterName rawValue varType i; 
settingsFile = 'configuration.mat';
save(settingsFile);
fprintf('Settings loaded from %s and saved in %s.\n', inputSettingsFile, settingsFile);
end