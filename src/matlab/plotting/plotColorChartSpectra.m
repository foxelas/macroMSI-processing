function [] = plotColorChartSpectra(x, vals, curveNames, options, fig)

     ylimits = [0, 2];
     hasMultiple = false;
     suffix2 = '';
     
    if ~isempty(options)
        hasMultiple = options{1};
    end 
    if length(options) > 1
        currentCase = options{2};
        parts = strsplit(currentCase, '_');
        currentCase = parts{1};
        if length(parts) > 1
            suffix2 = strcat(parts{2}, '_');
        end 
    end 

    if size(vals,2) == 36
        x = x; 
    elseif size(vals,2) == 17
        x = x(1,1:17);
    else
        x = x(1,18:end);
    end 
    
    suffix = '';
    if hasMultiple
        suffix = 'multipleAverage';
    end 
    figName = strcat(currentCase, 'calibrationSpectra');
    
    hasReflectanceRatio = true; 
    switch currentCase
        case 'expected'
            figTitle = 'Standard Spectra for Color Patches from Babel Color';
            ylimits  = [0,1];
        case 'measured'
            figTitle = 'Measured Spectra for Color Patches from our system';
        case 'measured-raw'
            figTitle = 'Raw Measured Spectra for Color Patches from our system';
            hasReflectanceRatio = false;
        case 'measured-adjusted'
            figTitle = 'Measured Spectra for Color Patches from our system after Adjustment';
            ylimits  = [0,1];
        case 'difference'
            figTitle = 'Percentage difference for Expected vs Measured spectra for Color Patches from our system';
            hasReflectanceRatio = false;
            ylimits = [-1,1];
        case 'difference-adjusted'
            figTitle = 'Percentage difference for Expected vs Measured spectra for Color Patches from our system after Adjustment';
            hasReflectanceRatio = false;
            ylimits = [-1,1];
        case 'points'
            figTitle = 'Spectra from random points of the image';
            hasReflectanceRatio = false;
            figName = strcat(getSetting('experiment') , '-pointSpectra');
        case 'measured-points'
            figTitle = 'Measured Spectra from our system';
        case '1x1window'
        case '2x2window'
        case '3x3window'
            parts = strsplit(currentCase, 'x');
            windowDim = str2num(parts(1));
            sprintf('Spectra from random points of the white image %dx%d spatial smoothing', windowDim, windowDim)
            hasReflectanceRatio = false;
            figName = strcat(getSetting('experiment') ,strcat('-pointSpectra', num2str(windowDim)));
        otherwise 
                error('Unsupported data name.');
    end 
    
    if contains(currentCase, 'difference')
        ylab = 'Difference of Reflectance Spectra (a.u.)';
    elseif hasReflectanceRatio 
        ylab = 'Reflectance Spectrum (%)';
    else 
        ylab = 'Reflectance Spectrum (a.u.)';
    end 
    
    if length(options) > 2
        ylimits = options{3};
    end 
    
    n = size(vals,1);
    if ~contains(currentCase, 'points')
        plotColors = getColorChartColors(curveNames);
    else 
        colors = hsv(n + 1);
        plotColors = colors(2:end, :);
    end
    
    
    if hasMultiple 
       %% Prepare Legend 
        h = zeros(n + 1,1);

        for i = 1:n 
            hold on;
            h(i) = plot(nan, nan, 'DisplayName', curveNames{i},'Color',  plotColors(i,:), 'LineWidth', 0.5);
            hold off;
        end 
        hold on;
        h(n+1) = plot(nan, nan, '--', 'DisplayName', 'Multiple Capture Average', 'Color', 'k', 'LineWidth', 1.5);
        hold off;
        
        %% Plot Spectra 
        for i = 1:n
            if hasMultiple 
                for k = 1:size(vals, 3)
                    hold on;
                    plot(x, squeeze(vals(i,:,k)), 'Color', plotColors(i,:), 'LineWidth', 0.5);
                    hold off;
                end 
                hold on;
                plot(x, squeeze(mean(vals(i,:,:), 3)), '--', 'Color', plotColors(i,:), 'LineWidth', 1.5);
                hold off;
            end 
        end
    else 
        h = zeros(n,1);
        for i = 1:n
            hold on
            h(i) = plot(x, vals(i,:), 'DisplayName', curveNames{i}, 'Color', plotColors(i,:), 'LineWidth', 3);
%             pause(0.1);
            hold off;
        end
    end
    
    xlim([420, 780]);
    ylim(ylimits);
    
    if hasReflectanceRatio 
        yline(1,'--','100%','LineWidth',3, 'DisplayName', 'Max Value');
    end
              
    xlabel('Wavelength (nm)', 'FontSize', 15);
    ylabel(ylab, 'FontSize', 15);
    title(figTitle, 'FontSize', 15);
    legend(h, 'Location', 'EastOutside');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

    plotName = fullfile(getSetting('savedir'), getSetting('saveFolder'), strcat(suffix2, figName, suffix, '.png'));
    setSetting('plotName', plotName);
    savePlot(fig);
end 