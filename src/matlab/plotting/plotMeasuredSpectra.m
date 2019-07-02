function [] = plotMeasuredSpectra(ID, Spectra, fig, saveOptions)

    if (nargin < 3)
        fig = figure;
    else 
        figure(fig);
        clf(fig);
    end

    if (nargin < 4)
        saveOptions.SaveImage = false;
    end 

    warning('off');
    wavelength = 380:5:780;
    average1 = [];
    average2 = [];
    n = size(Spectra, 1);
    for i = 1:n
    hold on
    if (ID(i).IsBenign)
        p1 = plot(wavelength, Spectra(i,:)*100, 'g', 'LineWidth', 3, 'DisplayName', 'Benign');
        average1 = [average1; Spectra(i,:)];
    else
        p2 = plot(wavelength, Spectra(i,:)*100, 'r', 'LineWidth', 3, 'DisplayName', 'Malignant');
        average2 = [average2; Spectra(i,:)];
    end
    hold off
    end
    
    legend([p1, p2], 'Location','northwest', 'FontSize', 15)
    xlabel('Wavelength (nm)', 'FontSize', 15);
    ylabel('Reflectance %', 'FontSize', 15);
    ylim([0,100]);
    saveOptions.saveImages = true;
    saveOptions.plotName = fullfile(saveOptions.savedir, 'general', 'allMeasuredSpectra');
    savePlot(fig, saveOptions);
        
    figure(fig+1);
    clf(fig+1);
    hold on
    e1 = errorbar(wavelength, mean(average1)*100, std(average1)/sqrt(n) *100,'-gs','MarkerSize',10,...
    'MarkerEdgeColor','green','MarkerFaceColor','green', 'DisplayName', 'Benign');
    e2 = errorbar(wavelength, mean(average2)*100, std(average2)/sqrt(n)*100,'-rs','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red', 'DisplayName', 'Malignant');
    e1.MarkerSize = 5;
    e1.Marker = 'o';
    e2.MarkerSize = 5;
    e2.Marker = 'o';

    hold off
    xlabel('Wavelength (nm)', 'FontSize', 15);
    ylabel('Reflectance %', 'FontSize', 15);
    ylim([0,40]);
    legend('FontSize', 15, 'Location','northwest')
    set(gcf, 'Position', get(0, 'Screensize'));
    saveOptions.plotName = fullfile(saveOptions.savedir, 'general', 'averageSpectra');
    savePlot(fig+1, saveOptions);
    
    warning('on');
end