function [spectralData, imageXYZ, wavelengths] = loadH5Data(filename)

saveFilename = fullfile(getSetting('matdir'),  strcat(filename, '.mat'));
if ~exist(saveFilename, 'file') 
    indir = getSetting('datadir');
    currentFile = fullfile(indir, filename); 
    %h5disp(currentFile);
    %h5info(currentFile);

    spectralData = h5read(currentFile, '/SpectralImage');
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