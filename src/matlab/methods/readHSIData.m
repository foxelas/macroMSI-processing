function [raw] = readHSIData(content, target, experiment, blackIsCapOn)
%%READHSIDATA returns the three images necessary for data analysis 
%   [raw] = readHSIData(target, experiment, blackIsCapOn)
 
if nargin < 4
    blackIsCapOn = false; 
end 

configuration = getSetting('configuration');
integrationTime = getSetting('integrationTime');
baseDir = fullfile(getSetting('matdir'), configuration, num2str(integrationTime));

%% Target image 
fcTarget = getFileConditions(content, target);
filename = getFilename(fcTarget{:});
[raw, ~, ~] = loadH5Data(filename, experiment);
figure(1); 
imshow(getDisplayImage(raw, 'rgb'));
setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('experiment'), num2str(getSetting('integrationTime')), target));
savePlot(1);
saveName = mkNewDir(strcat(baseDir, '_', target, '_target.mat')); 
save(saveName, 'raw', '-v7.3');

if ~strcmp(getSetting('normalization'), 'raw')
    %% White image 
    fcWhite = getFileConditions('whiteReflectance', target);
    filename = getFilename(fcWhite{:});
    [white, ~, wavelengths] = loadH5Data(filename, experiment);
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
    saveName = mkNewDir(strcat(baseDir, '_', target, '_white.mat')); 
    save(saveName, 'fullReflectanceByPixel', 'uniSpectrum', 'bandmaxSpectrum', '-v7.3');

    %% Black Image 
    if blackIsCapOn 
        fcBlack = getFileConditions('capOn', target);
    else 
        fcBlack = getFileConditions('lightsOff', target);
    end
    filename = getFilename(fcBlack{:});
    [blackReflectance, ~, ~] = loadH5Data(filename, experiment);
    figure(3); 
    imshow(getDisplayImage(blackReflectance, 'rgb'));
    setSetting('plotName', mkNewDir(getSetting('savedir'), getSetting('experiment'), strcat(target, '_black')));
    savePlot(3);

    saveName = mkNewDir(strcat(baseDir, '_', target,'_black.mat')); 
    save(saveName, 'blackReflectance', '-v7.3');
else 
    disp('Read only capture data, ignore white and black images.');
end 
end 