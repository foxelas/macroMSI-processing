function [spectralData, imageXYZ, wavelengths] = loadH5Data(fileName)

indir = getSetting('datadir');
currentFile = fullfile(indir, fileName);
h5disp(currentFile)
h5info(currentFile)

spectralData = h5read(currentFile, '/SpectralImage');
wavelengths = h5read(currentFile, '/Wavelengths');
imageX = h5read(currentFile, '/MeasurementImages/Tristimulus_X');
imageY = h5read(currentFile, '/MeasurementImages/Tristimulus_Y');
imageZ = h5read(currentFile, '/MeasurementImages/Tristimulus_Z');

imageXYZ = cat(3, imageX, imageY, imageZ);

save(strcat(fileName, '.mat'), 'spectralData', 'imageXYZ', 'wavelengths', '-v7.3');
end