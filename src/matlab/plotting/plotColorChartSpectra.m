function [] = plotColorChartSpectra(x, vals, spectraColorOrder, name, fig)

    
    ylab = 'Reflectance Spectrum (%)';
    ylimits = [0, 2];
    currentCase = name{2};
    switch currentCase
        case 'expected'
            figTitle = 'Standard Spectra for Color Patches from Babel Color';
            ylimits  = [0,1];
        case 'measured'
            figTitle = 'Measured Spectra for Color Patches from our system';
        case 'measured-raw'
            figTitle = 'Raw Measured Spectra for Color Patches from our system';
            ylab = 'Reflectance Spectrum (a.u.)';
        case 'measured-adjusted'
            figTitle = 'Measured Spectra for Color Patches from our system after Adjustment';
            ylimits  = [0,1];
        case 'difference'
            figTitle = 'Expected-Measured spectra for Color Patches from our system';
            ylab = 'Reflectance Difference (a.u,)';
        case 'difference-adjusted'
            figTitle = 'Expected-Measured spectra for Color Patches from our system after Adjustment';
            ylab = 'Reflectance Difference (a.u,)';
        otherwise 
            if numel(name) > 2
                figTitle = name{3};
                ylab = name{4};
                ylimits = name{5};
            else 
                error('Unsupported data name.');
            end 
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
        ylim(ylimits); %ylim([0, 1.2]);
        yline(1,'--','100%','LineWidth',3, 'DisplayName', 'Max Value');
    end
    
    if strcmp(ylab, 'Reflectance Difference (a.u,)')
        ylim([-1, 1]); 
    end
    
    xlabel('Wavelength (nm)', 'FontSize', 15);
    ylabel(ylab, 'FontSize', 15);
    title(figTitle, 'FontSize', 15);
    legend('Location', 'EastOutside');

    setSetting('plotName', fullfile(getSetting('savedir'), name{1}, strcat(currentCase, 'calibrationSpectra.png')));
    savePlot(fig);
end 