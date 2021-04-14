function [spectralData] = readHSIData(content, target, experiment, blackIsCapOn)
%%READHSIDATA returns the three images necessary for data analysis 
%   [raw] = readHSIData(target, experiment, blackIsCapOn)
 
if nargin < 4
    blackIsCapOn = false; 
end 

baseDir = mkNewDir(getSetting('matdir'), strcat(getSetting('database'), 'Triplets\'));

%% Target image 
fcTarget = getFileConditions(content, target);
[filename, tableId] = getFilename(fcTarget{:});
saveName = fullfile(baseDir, num2str(tableId));
[spectralData, ~, ~] = loadH5Data(filename);

if ~exist(strcat(saveName, '_target.mat'), 'file') || ~exist(strcat(saveName, '_black.mat'), 'file') || ~exist(strcat(saveName, '_white.mat'), 'file')
    figure(1); 
    imshow(getDisplayImage(spectralData, 'rgb'));
    setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('experiment'), num2str(getSetting('integrationTime')), target));
    savePlot(1);
    save(strcat(saveName, '_target.mat'), 'spectralData', '-v7.3');

    if ~strcmp(getSetting('normalization'), 'raw')
        %% White image 
        fcWhite = getFileConditions('whiteReflectance', target);
        filename = getFilename(fcWhite{:});
        [white, ~, wavelengths] = loadH5Data(filename);
        figure(2); 
        imshow(getDisplayImage(white, 'rgb'));
        setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('experiment'), strcat(target, '_white')));
        savePlot(2);

        %%UniSpectrum
        uniSpectrum = getSpectrumCurves(white);
        setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('experiment'), strcat(target, '_whitePlot_unispectrum')));
        plots(4, @plotSpectra, uniSpectrum, wavelengths, '99%-white', 'Reflectance Spectrum of White Balance Sheet');

        %%BandMax
        [m,n,w] = size(white);
        bandmaxSpectrum = max(reshape(white, m * n, w), [], 1);
        setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('experiment'), strcat(target, '_whitePlot_bandmax')));
        plots(5, @plotSpectra, bandmaxSpectrum, wavelengths, 'Bandmax spectrum', 'Bandmax Spectrum for the current Image');

        fullReflectanceByPixel = white;    
        save(strcat(saveName, '_white.mat'), 'fullReflectanceByPixel', 'uniSpectrum', 'bandmaxSpectrum', '-v7.3');

        %% Black Image 
        if blackIsCapOn 
            fcBlack = getFileConditions('capOn', target);
        else 
            fcBlack = getFileConditions('lightsOff', target);
        end
        filename = getFilename(fcBlack{:});
        [blackReflectance, ~, ~] = loadH5Data(filename);
        figure(3); 
        imshow(getDisplayImage(blackReflectance, 'rgb'));
        setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('experiment'), strcat(target, '_black')));
        savePlot(3);

        save(strcat(saveName, '_black.mat'), 'blackReflectance', '-v7.3');
    else 
        disp('Read only capture data, ignore white and black images.');
    end 
end 

end 