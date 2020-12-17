function [] = plotColorChartSpectra(x, vals, spectraColorOrder, name, fig)

    
    ylab = 'Reflectance Spectrum (%)';
    switch name 
        case 'expected'
            figTitle = 'Standard Spectra for Color Patches from Babel Color';
        case 'measured'
            figTitle = 'Measured Spectra for Color Patches from our system';
        case 'measured-raw'
            figTitle = 'Raw Measured Spectra for Color Patches from our system';
            ylab = 'Reflectance Spectrum (a.u.)';
        otherwise 
            error('Unsupported data name.');
    end 
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plotColors = getColorChartColors(spectraColorOrder);
    for i = 1: size(vals,1)
        hold on
        plot(x, vals(i,:), 'g', 'DisplayName', spectraColorOrder{i}, 'Color', plotColors(i,:), 'LineWidth', 3);
        pause(0.1);
        hold off;
    end
    
    if strcmp(ylab, 'Reflectance Spectrum (%)')
        xlim([420, 730]);
        ylim([0, 1.2]);
        yline(1,'--','100%','LineWidth',3, 'DisplayName', 'Max Value');
    end
    xlabel('Wavelength (nm)', 'FontSize', 15);
    ylabel(ylab, 'FontSize', 15);
    title(figTitle, 'FontSize', 15);
    legend('Location', 'EastOutside');

    setSetting('plotName', fullfile(getSetting('savedir'), strcat(name, '.png')));
    savePlot(fig);
end 