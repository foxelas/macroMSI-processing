function [spectralData, imageXYZ, wavelengths] = loadH5Data(filename, configuration)
%LOADH5DATA loads info from h5 file
%
%   [spectralData, imageXYZ, wavelengths] = loadH5Data(filename,
%   configuration) returns spectralData, XYZ image and capture wavelengths
%

saveFilename = mkNewDir(getSetting('matdir'), configuration, strcat(filename, '.mat'));

if ~exist(saveFilename, 'file')
    
    currentFile = adjustFilename(filename);
    %h5disp(currentFile);
    %h5info(currentFile);
    
    spectralData = gpuArray(h5read(currentFile, '/SpectralImage'));
    wavelengths = h5read(currentFile, '/Wavelengths');
    imageX = h5read(currentFile, '/MeasurementImages/Tristimulus_X');
    imageY = h5read(currentFile, '/MeasurementImages/Tristimulus_Y');
    imageZ = h5read(currentFile, '/MeasurementImages/Tristimulus_Z');
    imageXYZ = cat(3, imageX, imageY, imageZ);
    
    save(saveFilename, 'spectralData', 'imageXYZ', 'wavelengths', '-v7.3');
else
    load(saveFilename, 'spectralData', 'imageXYZ', 'wavelengths');
end

end

function currentFile = adjustFilename(filename)
indir = getSetting('datadir');
%     indir = 'C:\Users\elena\Desktop\'; %temporary read dir

filenameParts = strsplit(filename, '_');
dataDate = filenameParts{1};
if ~contains(indir, dataDate)
    filenameParts = strsplit(indir, '\\saitama');
    ending = filenameParts{2};
    filenameParts = strsplit(ending, '_');
    oldDate = filenameParts{1};
    indir = strrep(indir, oldDate, dataDate);
end

currentFile = fullfile(indir, filename);
end