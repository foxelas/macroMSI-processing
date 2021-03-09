
%% Initialization
close all;

userSettingsFile = '..\..\conf\hsiUserSettings.csv';

originDir = 'F:\temp\mspi';
configuration = 'singleLightClose';

%% Main
setOpt(userSettingsFile);

%% Other settings
if ~exist('dataDate', 'var')
    dataDate = '20210127';
    warning('Setting default date: 20210127.');
end

indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'), 'h5');
if exist('indirFolder', 'var') && ~isempty(indirFolder)
    indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'), indirFolder, 'h5');
end

setSetting('datadir', indir);
matdir = fullfile(originDir, 'matfiles\hsi');
setSetting('matdir', matdir);

if exist('experiment', 'var')
    setSetting('experiment', experiment);
    setSetting('saveFolder', experiment);
end

if exist('integrationTime', 'var')
    setSetting('integrationTime', integrationTime);
end

if ~exist('configuration', 'var')
    configuration = 'singleLightClose';
end
setSetting('configuration', configuration);

if ~exist('normByPixel', 'var')
    normByPixel = true;
end
setSetting('normByPixel', normByPixel);

if exist('targetPosition', 'var')
    setSetting('targetPosition', targetPosition);
end