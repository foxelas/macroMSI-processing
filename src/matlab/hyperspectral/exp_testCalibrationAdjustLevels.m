%% Compare values on the white image and the error between them
% compare with the measured values from source closer to the light source
% and calculate integration coefficients to adjust values to same level 

%% Setup 
startRun;

%% Colorchart with different normalizations and positions
experiment = 'testCalibrationAdjustLevels';
setSetting('experiment', experiment);
dataDate = '20210308';
configuration = 'singleLightClose';
integrationTime = 1460;
normalization = 'raw';
version = '1';
initialization;

%% Read White image as base 
setSetting('saveFolder', experiment);
savedir = fullfile(getSetting('savedir'), getSetting('saveFolder'));
target = strcat('stomachTissue', version);
fileConditions = getFileConditions('whiteReflectance', target);
rawWhite = readHSIData('whiteReflectance', target, experiment);
[m,n,w] = size(rawWhite);
xPoints = [100, floor(m/2) + 100, (m - 100)];
yPoints = [100, floor(n/2) + 100, (n - 100)];
[measured, curveNames] = getRepresentativePoints(fileConditions, xPoints, yPoints, [0, 0.003]);

rawTissue = readHSIData('tissue', target, experiment);

%% Find image point most responsible for bandmax 
wavelengths = getWavelengths(w);
[bandmaxSpectrum, maxIds] = max(rawWhite, [], [1,2], 'linear');
%ignore indexes below 420nm because of noise
maxIds = maxIds(40:end);
[mX, mY,~] = arrayfun(@(x) ind2subArray(m,n,w,x), squeeze(maxIds));
setSetting('plotName', mkNewDir(savedir, strcat(target, '_whitePlot_bandmax')));
plots(5, @plotSpectra, bandmaxSpectrum, wavelengths, 'Bandmax spectrum', 'Bandmax Spectrum for the current Image');

a = unique(maxIds);
hc = [a,histc(maxIds(:),a)];
freqInd = find(hc(:, 2) == max(hc(:,2)));
[row, col] = ind2sub([m, n], a(freqInd));
maxfreq = hc(freqInd, 2);
if length(maxfreq) == 1
    fprintf('Pixel (%d, %d) has bandmax value %d times (most common).\n', col, row, maxfreq);
else 
    row = floor(mean(mX));
    col = floor(mean(mY));
    fprintf('Average position (%d, %d) of pixels contributing montsly to bandmax spectral values.\n', col, row);
end

dispImage = getDisplayImage(rawWhite);
fig = figure(6);clf;
imshow(dispImage);
hold on;
for i = 1:length(mY)
    plot(mX(i), mY(i), 'bd', 'MarkerSize', 10, 'LineWidth', 5);
end
plot(col, row, 'rd', 'MarkerSize', 10, 'LineWidth', 5);
hold off;
setSetting('plotName', mkNewDir(savedir, strcat(target, '_bandmax_participants')));
savePlot(fig);
%%Now this pixel will be used as reference and is spectrum as reference spectrum
%% Linear Regression
yVals = squeeze(rawWhite(row, col, :));
xVals = reshape(rawWhite, [m*n,w]);
lvlAdjCoeff =   repmat(yVals', [m*n,1]) ./ xVals;

lvlAdjCoeff = reshape(lvlAdjCoeff, [m,n,w]);

fig = figure(1);clf;
b = permute(lvlAdjCoeff(:,:,100), [2, 1]);
imagesc(b);
colorbar;
title('At 480nm');
setSetting('plotName', fullfile( savedir, 'lr_coeffs_480.png'));
savePlot(fig);

fig = figure(2);clf;
b =  permute(lvlAdjCoeff(:,:,350), [2, 1]);
imagesc(b);
colorbar;
title('At 730nm');
setSetting('plotName', fullfile( savedir, 'lr_coeffs_730.png'));
savePlot(fig);

fig = figure(3);clf;
v3 = permute(mean(lvlAdjCoeff, 3), [2,1]);
imagesc(v3)
colorbar;
title('Average multiplier to adjust to center of screen level');
setSetting('plotName', fullfile( savedir, 'lr_coeffs_average.png'));
savePlot(fig);

%% Now plot the adjusted values for white reflectance 
setSetting('saveFolder', fullfile(getSetting('experiment'), 'lvlAdj'));
[measured, curveNames] = getRepresentativePoints(fileConditions, xPoints, yPoints, [0, 0.003], [], lvlAdjCoeff);

%% Now read the original normalized spectra (for white)
setSetting('saveFolder', experiment);
setSetting('normalization', 'byPixel');
setSetting('saveFolder', getSetting('experiment'));
fileConditions = getFileConditions('whiteReflectance', target);
targetName = fileConditions{4};
spectralData = normalizeHSI(rawWhite, targetName);

[m,n,w] = size(rawWhite);
xPoints = [100, floor(m/2) + 100, (m - 100)];
yPoints = [100, floor(n/2) + 100, (n - 100)];
setSetting('plotName', fullfile(savedir, strcat(target,'_', 'white')));
[measured, curveNames] = getRepresentativePoints(spectralData, xPoints, yPoints, [0, 1.2]);

%% After applying the result  (for white)
setSetting('saveFolder', fullfile(getSetting('experiment'), 'lvlAdj'));
setSetting('plotName', fullfile(savedir, strcat(target,'_', 'white')));
[measured, curveNames] = getRepresentativePoints(spectralData, xPoints, yPoints, [0, 1.2], [], lvlAdjCoeff);


%% Now read the original normalized spectra 
setSetting('saveFolder', getSetting('experiment'));

fileConditions = getFileConditions('tissue', target);
targetName = fileConditions{4};
spectralData = normalizeHSI(rawTissue, targetName);

xPoints = [250, 320, 400];
yPoints = [250, 350, 400];
setSetting('plotName', fullfile(savedir, strcat(target,'_', 'tissue')));
[measured, curveNames] = getRepresentativePoints(spectralData, xPoints, yPoints, [0, 1.2]);

%% After applying the result 
setSetting('saveFolder', fullfile(getSetting('experiment'), 'lvlAdj'));
setSetting('plotName', fullfile(savedir, strcat(target,'_', 'tissue')));
[measured, curveNames] = getRepresentativePoints(spectralData, xPoints, yPoints, [0, 1.2], [], lvlAdjCoeff);

endRun;


function [ii,jj,kk]  = ind2subArray(a,b,c,y)
[ii,jj,kk] = ind2sub([a,b,c], y);
end