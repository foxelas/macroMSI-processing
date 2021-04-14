
%% Initialization
disp('Initialization started');

userSettingsFile = '..\..\conf\hsiUserSettings.csv';
originDir = 'D:\elena\mspi';

%% Main
setOpt(userSettingsFile);

%% Other settings
if ~exist('dataDate', 'var')
    dataDate = '20210127';
    warning('Setting default date: 20210127.');
end
setSetting('dataDate', dataDate);

indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'), 'h5');
if exist('indirFolder', 'var') && ~isempty(indirFolder)
    indir = fullfile(originDir, '2_saitamaHSI', strcat('saitama', dataDate, '_test'), indirFolder, 'h5');
end

if ~exist('database', 'var')
    dataDate = 'calib';
    warning('Setting default database for calibration: calib.');
end
setSetting('database', database);

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

if exist('normalization', 'var')
    setSetting('normalization', normalization);
end

disp('Initialization finished');

