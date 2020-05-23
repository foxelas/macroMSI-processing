function [] = plotMeasuredSpectra(ID, Spectra, fig)

warning('off');

%% All measured
wavelength = 380:5:780;
average1 = [];
average2 = [];
averageFixed = [];
averageUnfixed = [];
n = size(Spectra, 1);
for i = 1:n
    hold on
    if (ID(i).IsBenign)
        p1 = plot(wavelength, Spectra(i, :)*100, 'g', 'LineWidth', 3, 'DisplayName', 'Benign');
        average1 = [average1; Spectra(i, :)];
    else
        p2 = plot(wavelength, Spectra(i, :)*100, 'r', 'LineWidth', 3, 'DisplayName', 'Malignant');
        average2 = [average2; Spectra(i, :)];
    end
    hold off
    if (ID(i).IsFixed)
        averageFixed = [averageFixed; Spectra(i, :)];
    else
        averageUnfixed = [averageUnfixed; Spectra(i, :)];
    end

end

legend([p1, p2], 'Location', 'northwest', 'FontSize', 15)
xlabel('Wavelength (nm)', 'FontSize', 15);
ylabel('Reflectance %', 'FontSize', 15);
ylim([0, 100]);
setSetting('plotName', fullfile(getSetting('savedir'), 'general', 'allMeasuredSpectra'));
savePlot(fig);

%% Malignancy average
figure(fig+1);
clf(fig+1);
hold on
e1 = errorbar(wavelength, mean(average1)*100, std(average1)/sqrt(n)*100, '-gs', 'MarkerSize', 10, ...
    'MarkerEdgeColor', 'green', 'MarkerFaceColor', 'green', 'DisplayName', 'Benign');
e2 = errorbar(wavelength, mean(average2)*100, std(average2)/sqrt(n)*100, '-rs', 'MarkerSize', 10, ...
    'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red', 'DisplayName', 'Malignant');
e1.MarkerSize = 5;
e1.Marker = 'o';
e2.MarkerSize = 5;
e2.Marker = 'o';
hold off
xlabel('Wavelength (nm)', 'FontSize', 20);
ylabel('Reflectance %', 'FontSize', 20);
ylim([0, 40]);
legend('FontSize', 20, 'Location', 'northwest')
set(gcf, 'Position', get(0, 'Screensize'));
setSetting('plotName', fullfile(getSetting('savedir'), 'general', 'averageSpectra'));
savePlot(fig+1);

%% Fixing average
figure(fig+1);
clf(fig+1);
hold on
e1 = errorbar(wavelength, mean(averageFixed)*100, std(average1)/sqrt(n)*100, '-gs', 'MarkerSize', 10, ...
    'MarkerEdgeColor', 'green', 'MarkerFaceColor', 'green', 'DisplayName', 'Fixed');
e2 = errorbar(wavelength, mean(averageUnfixed)*100, std(average2)/sqrt(n)*100, '-rs', 'MarkerSize', 10, ...
    'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red', 'DisplayName', 'Unfixed');
e1.MarkerSize = 5;
e1.Marker = 'o';
e2.MarkerSize = 5;
e2.Marker = 'o';

hold off
xlabel('Wavelength (nm)', 'FontSize', 20);
ylabel('Reflectance %', 'FontSize', 20);
ylim([0, 40]);
legend('FontSize', 20, 'Location', 'northwest')
set(gcf, 'Position', get(0, 'Screensize'));
setSetting('plotName', fullfile(getSetting('savedir'), 'general', 'averageSpectraFixing'));
savePlot(fig+1);

benId = 85;
setSetting('plotName', fullfile(getSetting('savedir'), 'general', 'benunf.png'));
plotFunWrapper(1, @plotReconstructedCurves, Spectra(benId, :)', {'Benign, Unfixed'}, 380:5:780, '');
malId = 84;
setSetting('plotName', fullfile(getSetting('savedir'), 'general', 'malunf.png'));
plotFunWrapper(1, @plotReconstructedCurves, Spectra(malId, :)', {'Malignant, Unfixed'}, 380:5:780, '');
fixedBenId = 79;
setSetting('plotName', fullfile(getSetting('savedir'), 'general', 'benf.png'));
plotFunWrapper(1, @plotReconstructedCurves, Spectra(fixedBenId, :)', {'Benign, Fixed'}, 380:5:780, '');
fixedMalId = 78;
setSetting('plotName', fullfile(getSetting('savedir'), 'general', 'malf.png'));
plotFunWrapper(1, @plotReconstructedCurves, Spectra(fixedMalId, :)', {'Malignant, Fixed'}, 380:5:780, '');

warning('on');
end