function [] = plotReconstructedCurves(reflectanceSpectra, curveNames, wavelength, figTitle, ...
    markers, fig)

%{
		% returns a comparative plot of the estimation reflectanceSpectra versus the measured spectrum.
		%wavelength: vector, e.g. 401x1
		%reflectanceSpectra: the estimation results in order { 'Measured', 'Est-MSgreen', 'Est-MSrms','Est-MSadjusted', 'Est-MSextended', 'Est-RGB'}
		%name: the name of the sample, is included in the plot title
		%method: the estimation method (eg SmoothingMatrixMethod)
		%plotOutName = if it's not empty, then the plot is saved with that filename
		% msiReflectances = presents also a subplot os the MSI reflectances value
    %}

    fc = [450, 465, 505, 525, 575, 605, 630];


    % each column of 'spectrum' is that data for a plot line
    curveN = size(reflectanceSpectra, 2);

    if isempty(curveNames)
        curveNames = {'MS center \lambda', 'Measured', 'Est-MSgreen', 'Est-MSrms', 'Est-MSadjusted', 'Est-MSextended', 'Est-RGB'};
    else
        curveNames = cellfun(@(x) strrep(x, '_', ' '), curveNames, 'un', 0);
    end

    if ~exist('markers', 'var') || isempty(markers)
        for i = 1:curveN
            markers{i} = 'none';
        end
    end

    color = colorcube(curveN+3);
    color = color(3:end-1, :);

    [~, peakIdx] = ismember(fc, wavelength); % mark filter wavelengths

    hold on
    for i = 1:curveN
        currentCurvName = curveNames{i};
        if contains(lower(currentCurvName), 'rgb')
            lineStyle = ':';
            clrs = 'k';
        elseif contains(lower(currentCurvName), 'measured')
            lineStyle = '--';
            clrs = 'g';
        else
            lineStyle = '-';
            clrs = color(i, :);
        end
        plot(wavelength, reflectanceSpectra(:, i).*100, 'Color', clrs, 'Marker', markers{i}, ...
            'LineWidth', 2, 'LineStyle', lineStyle, 'DisplayName', currentCurvName); % plot estimated reflectances
    end
    plot(wavelength(peakIdx), reflectanceSpectra(peakIdx, 1).*100, 'rx', 'DisplayName', 'Channel Wavelength', 'LineWidth', 1.3); % plot measured reflectance
    hold off

    xlabel('Wavelength \lambda (nm)', 'FontSize', 15);
    ylabel('Reflectance %', 'FontSize', 15);
    xlim([400, 700]);
    %figTitle = strjoin({'Comparative plot of Wiener estimation for Sample', simpleSampleName(figTitle)}, ' ');
    %title(figTitle);
    legend([curveNames, 'Channel Wavelength'], 'Location', 'EastOutside', 'FontSize', 15); % 'Orientation','horizontal');

    %Optional
    set(gcf, 'Position', get(0, 'Screensize'));

    setSetting('cropBorders', true);
    savePlot(fig);

end

function [sampleName] = simpleSampleName(name)
attr = strsplit(name);
sampleName = attr{1};
end