%% Prepare rgb preview images for each  hsm file in the db 

%% Setup 
startRun;

%% Read 
dataDate = '20210317';
integrationTime = 1360;
configuration = 'singleLightClose';
normalization = 'byPixel';

initialization;

basedirParts = strsplit(getSetting('datadir'), '\');
basedirParts = basedirParts(1:(end-2));
basedir = fullfile(basedirParts{:});

dataTable = readtable(fullfile(getSetting('datasetSettingsDir'), 'testDB.xlsx'));
fig = figure(1);

for i=9:height(dataTable)
    currentRow = dataTable(i,:);
    filename = char(currentRow.Filename);
    folder = strcat('saitama',num2str(currentRow.CaptureDate), '_test');
    imgName = fullfile(basedir, folder, 'h5', strcat(strrep(filename, '.', '_'), '_preview.jpg'));
    
    if ~exist(imgName, 'file')
        [raw, ~, ~] = loadH5Data(filename, getSetting('experiment'));
        z = size(raw, 3) ;
        if z < 401 
            baseImage = rescale(squeeze(raw(:,:,min(z,100)))); 
        else 
            baseImage = getDisplayImage(raw, 'rgb');
        end
        setSetting('plotName', imgName);
        clf(fig);
        imshow(baseImage);
        title(strrep(filename, '_', '-'));
        savePlot(fig);
    end 
end 
