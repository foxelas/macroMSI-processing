webdir = getSetting('webdir');
eumelaninFilename = 'eumelanin_absroption.csv';
hbFilename = 'hb_absorption_spectra_prahl.csv';
eumelaninData = delimread(fullfile(webdir, eumelaninFilename), ',', 'num');
eumelaninData = eumelaninData.num;
hbData = delimread(fullfile(webdir, hbFilename), ',', 'num');
hbData = hbData.num;

eumelaninLambda = eumelaninData(:,1);
eumelanin1 = eumelaninData(:, 2);
eumelanin2 = eumelaninData(:, 3);

hbLambda = hbData(:,1);
hbO2 = hbData(:, 2);
hb = hbData(:, 3);

fig = figure(1);
clf;
hold on; 
%plot(eumelaninLambda, log10(eumelanin1), 'DisplayName', 'Eumelanin1', 'LineWidth', 2); %cm-1 / (mg/ml)
plot(eumelaninLambda, log10(eumelanin2), 'DisplayName', 'Eumelanin2', 'LineWidth', 2); %cm-1 / (moles/liter)
plot(hbLambda, log10(hbO2), 'DisplayName', 'HbO2', 'LineWidth', 2); %cm-1/M
plot(hbLambda, log10(hb), 'DisplayName', 'Hb', 'LineWidth', 2); %cm-1/M
hold off 
xlabel('Wavelength (nm)', 'FontSize', 15);
ylabel('log_{10}(Absorption) (cm^{-1}/ M)', 'FontSize', 15);
l = legend('Location', 'northeastoutside');
l.FontSize = 13;

hold on;
for i = 2:length(bands)
    xline(bands(i),'--', 'DisplayName', strcat('LED@', num2str(bands(i)), 'nm') );
end
hold off; 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);

setSetting('saveInHQ', true);
setSetting('plotName', fullfile(getSetting('savedir'), getSetting('common'), 'skinChromophoreAbsorption'));
savePlot(fig);



