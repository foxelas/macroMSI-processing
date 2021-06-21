if false 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for April 6th %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup 
startRun;

%% Read tissue 2
version = '3';
experiment = strcat('testStomach', version, 'embc');
dataDate = '20210406';
normalization = 'byPixel';
initialization;

targets = { strcat('stomachTissue', version), strcat('stomachTissue', version, '_int_auto_gain_4')};
types = {'tissue', 'tissue'};


target = 'fullscreen';
setSetting('normalization', 'raw');
integrationTimes = 618; %[1380, 800, 2130, 618];
m = length(integrationTimes);
measured = cell(m,1);
xPoints = [100, 650, 1200];
yPoints = [100, 500, 900];
for i = 1:m
    setSetting('integrationTime', integrationTimes(i));
    setSetting('saveFolder', 'comparsionOfWhiteRawMeasumerements');
    fileConditions = getFileConditions('whiteReflectance', target);
    [measured{i}, curveNames] = getRepresentativePoints(fileConditions, xPoints, yPoints);
end 
end 

%D:\elena\Google Drive\titech\research\experiments\output\2. Short Papers\4. Calibration\3_main\img 

if false 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for march 17 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
version = '2';
experiment = strcat('testStomach', version, 'embc2');
dataDate = '20210317';
normalization = 'byPixel';
initialization;

targets = { strcat('stomachTissue', version), strcat('stomachTissue', version, '_int_auto_gain_4')};
types = {'tissue', 'tissue'};


target = 'fullscreen';
setSetting('normalization', 'raw');
integrationTimes = 1360; %[1380, 800, 2130, 618];
m = length(integrationTimes);
measured = cell(m,1);
xPoints = [100, 650, 1200];
yPoints = [100, 500, 900];
for i = 1:m
    setSetting('integrationTime', integrationTimes(i));
    setSetting('saveFolder', 'comparsionOfWhiteRawMeasumerements');
    fileConditions = getFileConditions('whiteReflectance', target);
    [measured{i}, curveNames] = getRepresentativePoints(fileConditions, xPoints, yPoints);
end 
end 

if false 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% for april 4th tissue %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

version = '3';
experiment = strcat('testStomach', version);
dataDate = '20210406';
initialization;

normalizations = {'bandmax', 'uniSpectrum', 'byPixel', 'raw'};
m = numel(normalizations);
xPoints = [500, 600];
yPoints = [500];

for k = 1:m
        target = strcat('stomachTissue', version);
        setSetting('integrationTime', 618)
        normalization = normalizations{k};
        setSetting('saveFolder', 'comparsionOfNormalizationsMeasumerements');
        setSetting('normalization', normalization);    
        fileConditions = getFileConditions('tissue', target);
        [meas, curveNames] = getRepresentativePoints(fileConditions, xPoints, yPoints);
        measured(k,1:2,1:401) = squeeze(meas);
end 
v = fullfile(getSetting('savedir'),getSetting('saveFolder'), 'spectra.mat');
save(v, 'measured', 'normalizations', 'xPoints', 'yPoints', 'dataDate', 'experiment');

close all; 
%%%%%%%%%%%%%%%%%%%%
normalizations = {'Bandmax', 'Average', 'Pixelwise', 'raw'};
fig = figure(1);
x = getWavelengths(401);
plotColors = {'m', 'c', 'b'};
h = zeros(6,1);
% h = zeros(3,1);
hold on;
for i = 1:3 
    spectrum = squeeze(measured(i, 1, :));
    [spectrum, x] = CutRange(spectrum, x);
    h(i) = plot(x, spectrum, '-', 'DisplayName', strcat(normalizations{i}, '@P1'), 'Color', plotColors{i}, 'LineWidth', 2);    
%     h(i) = plot(x, spectrum, '-', 'DisplayName', normalizations{i}, 'Color', plotColors{i}, 'LineWidth', 2);    
end 
for i = 1:3 
    spectrum = squeeze(measured(i, 2, :));
    [spectrum, x] = CutRange(spectrum, x);
    h(i+3) = plot(x, spectrum, ':', 'DisplayName', strcat(normalizations{i}, '@P2'), 'Color', plotColors{i}, 'LineWidth', 1);    
end 
hold off; 
 
xlim([400, 750]);
ylim([0,1]);
xlabel('Wavelength (nm)', 'FontSize', 15);
ylab = 'Reflectance Spectra (a.u.)';
ylabel(ylab, 'FontSize', 15);
legend(h, 'Location', 'NorthWest');
ax = gca;
ax.YAxis.Exponent = 0;

plotName = fullfile(getSetting('savedir'), getSetting('saveFolder'), strcat('normalizations', '.png'));
setSetting('plotName', plotName);
savePlot(fig); 
 
%%%%%%%%%%%%%%%%%%%%
fig = figure(2);
h = zeros(2,1);
colors = {'c', 'm'};hold on 
for i = 1:2 
spectrum = squeeze(measured(4, i, :));
[spectrum, x] = CutRange(spectrum, x);
h(i) = plot(x, spectrum, '-', 'DisplayName', strcat('Raw@P', num2str(i)), 'Color', colors{i}, 'LineWidth', 2);    
end
hold off 
xlim([400, 750]);
ylim([0,0.0015]);
xlabel('Wavelength (nm)', 'FontSize', 15);
ylab = 'Measured Spectra (a.u.)';
ylabel(ylab, 'FontSize', 15);
legend('Location', 'NorthWest');
ax = gca;
ax.YAxis.Exponent = 0;
yticks(0:0.0005:0.0015)

plotName = fullfile(getSetting('savedir'), getSetting('saveFolder'), strcat('raw', '.png'));
setSetting('plotName', plotName);
savePlot(fig); 

%%%%%%%%%%%%%%%%%%%%
alpha = 0.595324592687991; 
fig = figure(3);
spectrum = squeeze(measured(3, 1, :)) * alpha;
[spectrum, x] = CutRange(spectrum, x);
plot(x, spectrum, '-', 'DisplayName', 'Adjusted', 'Color', 'c', 'LineWidth', 2);    
xlim([400, 750]);
ylim([0,1]);
xlabel('Wavelength (nm)', 'FontSize', 15);
ylab = 'Reflectance Spectra (a.u.)';
ylabel(ylab, 'FontSize', 15);
legend('Location', 'NorthWest');
ax = gca;
ax.YAxis.Exponent = 0;

plotName = fullfile(getSetting('savedir'), getSetting('saveFolder'), strcat('adjusted', '.png'));
setSetting('plotName', plotName);
savePlot(fig); 

%%%%%%%%%%%%%%%%%%%%
fig = figure(4);
h = zeros(2,1);
colors = {'c', 'm'};
hold on 
for i = 1:2 
spectrum = squeeze(measured(3, i, :)) * alpha;
[spectrum, x] = CutRange(spectrum, x);
plot(x, spectrum, '-', 'DisplayName', strcat('Adjusted@P', num2str(i)), 'Color', colors{i}, 'LineWidth', 2);    
end
hold off 
xlim([400, 750]);
ylim([0,1]);
xlabel('Wavelength (nm)', 'FontSize', 15);
ylab = 'Reflectance Spectra (a.u.)';
ylabel(ylab, 'FontSize', 15);
legend('Location', 'NorthWest');
ax = gca;
ax.YAxis.Exponent = 0;

plotName = fullfile(getSetting('savedir'), getSetting('saveFolder'), strcat('adjusted2', '.png'));
setSetting('plotName', plotName);
savePlot(fig); 

% plots(1, @plotColorChartSpectra, squeeze(measured(:, 1, :)), normalizations, 'measured', ...
%     [0, 1], false);
% 
% plots(2, @plotColorChartSpectra, squeeze(measured(:, 2, :)), normalizations, 'measured', ...
%     [0, 1], false);
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% new %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if true 
%     basedir = 'D:\elena\mspi\2_saitamaHSI\saitama20210127_test\h5\'; 
%     basedir2 = 'D:\elena\mspi\2_saitamaHSI\saitama20210107_test\h5\'; 
%     spectralDataFilter = readHSItemp(fullfile(basedir, '20210127_174840_wet_filter_1460.h5'));
%     white = readHSItemp(fullfile(basedir, '20210127_181706_white_filter_800.h5'));
%     black= readHSItemp(fullfile(basedir2, '20210107_211617_black.h5'));
%     baseImage = getDisplayImage(spectralDataFilter, 'rgb');
    xPoints = [477, 591, 370]; 
    yPoints = [636, 738, 401];
    plots(1, @plotPointsOnImage, baseImage, xPoints, yPoints, false);
    
    x = yPoints(1);
    y = xPoints(1);
    norm1 = squeeze((spectralDataFilter(x,y,:) - black(x,y,:)) ./ (white(x,y,:) - black(x,y,:)));
    x = yPoints(2);
    y = xPoints(2);
    norm2 = squeeze((spectralDataFilter(x,y,:) - black(x,y,:)) ./ (white(x,y,:) - black(x,y,:))); 
    x = yPoints(3);
    y = xPoints(3);
    norm3 = squeeze((spectralDataFilter(x,y,:) - black(x,y,:)) ./ (white(x,y,:) - black(x,y,:))); 
    figure(2)
    w = getWavelengths(401);
    hold on 
    [norm1, w1] = CutRange(norm1, w);
    plot(w1, norm1);
    [norm2, ~] = CutRange(norm2, w);
    plot(w1, norm2);
    [norm3, ~] = CutRange(norm3, w);
    plot(w1, norm3);
    hold off 

    v = fullfile('D:\elena\Google Drive\titech\research\experiments\output\hsi\comparisonOfFilter', 'spectraFilter.mat');
    save(v, 'w1', 'norm1', 'norm2', 'norm3', 'xPoints', 'yPoints');

%     spectralDataNoFilter = readHSItemp(fullfile(basedir, '20210127_174840_wet_nofilter_300.h5'));
%     whiteNoFilter = readHSItemp(fullfile(basedir, '20210127_181706_white_nofilter_300.h5'));
%     blackNoFilter = readHSItemp(fullfile(basedir2, '20210107_211325_black.h5'));
%     baseImage2 = getDisplayImage(spectralDataNoFilter, 'rgb');
    xPoints = [562, 684, 492]; 
    yPoints = [674, 770, 480];
    plots(3, @plotPointsOnImage, baseImage2, xPoints, yPoints, false);
    
    x = yPoints(1);
    y = xPoints(1);
    norm1b = squeeze((spectralDataFilter(x,y,:) - black(x,y,:)) ./ (white(x,y,:) - black(x,y,:)));
    x = yPoints(2);
    y = xPoints(2);
    norm2b = squeeze((spectralDataFilter(x,y,:) - black(x,y,:)) ./ (white(x,y,:) - black(x,y,:))); 
    x = yPoints(3);
    y = xPoints(3);
    norm3b = squeeze((spectralDataFilter(x,y,:) - black(x,y,:)) ./ (white(x,y,:) - black(x,y,:))); 
    figure(4)
    w = getWavelengths(401);
    hold on 
    [norm1b, w1] = CutRange(norm1b, w);
    plot(w1, norm1b);
    [norm2b, ~] = CutRange(norm2b, w);
    plot(w1, norm2b);
    [norm3b, ~] = CutRange(norm3b, w);
    plot(w1, norm3b);
    hold off 

    v = fullfile('D:\elena\Google Drive\titech\research\experiments\output\hsi\comparisonOfFilter', 'spectraNoFilter.mat');
    save(v, 'w1', 'norm1b', 'norm2b', 'norm3b', 'xPoints', 'yPoints');

    figure(5);
    hold on 
    plot(w1, norm1, '-', 'DisplayName', 'With Filter @P1', 'Color', 'm', 'LineWidth', 2);
    plot(w1, norm2, '-', 'DisplayName', 'With Filter @P2', 'Color', 'c', 'LineWidth', 2);
    plot(w1, norm3, '-', 'DisplayName', 'With Filter @P3', 'Color', 'b', 'LineWidth', 2);
    plot(w1, norm1b, '--', 'DisplayName', 'No Filter @P1', 'Color', 'm', 'LineWidth', 1);
    plot(w1, norm2b, '--', 'DisplayName', 'No Filter @P2', 'Color', 'c', 'LineWidth', 1);
    plot(w1, norm3b, '--', 'DisplayName', 'No Filter @P3', 'Color', 'b', 'LineWidth', 1);
    hold off 
    xlim([400, 750]);
    ylim([0,1.5]);
    xlabel('Wavelength (nm)', 'FontSize', 15);
    ylab = 'Reflectance Spectra (a.u.)';
    ylabel(ylab, 'FontSize', 15);
    legend('Location', 'NorthWest');
    ax = gca;
    ax.YAxis.Exponent = 0;

    figure(6);
    hold on 
    plot(w1, norm1, '-', 'DisplayName', 'With Filter @P1', 'Color', 'm', 'LineWidth', 2);
    plot(w1, norm2, '-', 'DisplayName', 'With Filter @P2', 'Color', 'c', 'LineWidth', 2);
    plot(w1, norm3, '-', 'DisplayName', 'With Filter @P3', 'Color', 'b', 'LineWidth', 2);
    hold off 
    xlim([400, 750]);
    ylim([0,1.5]);
    xlabel('Wavelength (nm)', 'FontSize', 15);
    ylab = 'Reflectance Spectra (a.u.)';
    ylabel(ylab, 'FontSize', 15);
    legend('Location', 'NorthWest');
    ax = gca;
    ax.YAxis.Exponent = 0;

    figure(7);
    hold on 
    plot(w1, norm1b, '-', 'DisplayName', 'No Filter @P1', 'Color', 'm', 'LineWidth', 2);
    plot(w1, norm2b, '-', 'DisplayName', 'No Filter @P2', 'Color', 'c', 'LineWidth', 2);
    plot(w1, norm3b, '-', 'DisplayName', 'No Filter @P3', 'Color', 'b', 'LineWidth', 2);
    hold off 
    xlim([400, 750]);
    ylim([0,1.5]);
    xlabel('Wavelength (nm)', 'FontSize', 15);
    ylab = 'Reflectance Spectra (a.u.)';
    ylabel(ylab, 'FontSize', 15);
    legend('Location', 'NorthWest');
    ax = gca;
    ax.YAxis.Exponent = 0;
end 

function spectralData = readHSItemp(currentFile)
spectralData = double(h5read(currentFile, '/SpectralImage'));
end 