function [] = plotReconstructedCurves(reflectanceSpectra, curveNames, wavelength,  figTitle, ...
      fig,saveOptions, markers)

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
    
	if (nargin < 5)
		fig = figure;
	else     
		figure(fig);
        clf(fig);
	end
	if (nargin < 6)
		saveOptions.SaveImage = false;
	end

	% each column of 'spectrum' is that data for a plot line
	curveN = size(reflectanceSpectra, 2);
    
	if isempty(curveNames)
		curveNames = {'MS center \lambda', 'Measured', 'Est-MSgreen', 'Est-MSrms', 'Est-MSadjusted', 'Est-MSextended', 'Est-RGB'};
	else
		curveNames = cellfun(@(x) strrep(x, '_', ' '), curveNames, 'un', 0);
	end
	
	if ~exist('markers', 'var') ||  isempty(markers)
		for i = 1:curveN
			markers{i} = 'none';
		end
	end
		
	color = colorcube(curveN+10);         
	[~, peakIdx] = ismember(fc, wavelength); % mark filter wavelengths
	
	hold on
	for i = 1:curveN
        if (curveN > 1); currentCurvName = curveNames{i + 1}; else; currentCurvName= curveNames{i}; end
		if contains( lower(currentCurvName), 'rgb')
			lineStyle = ':';
        elseif contains( lower(currentCurvName), 'measured')
			lineStyle = '--';
		else 
			lineStyle = '-';
        end
        if (curveN > 1); clrs = color(i, :); else; clrs = 'b'; end
		plot(wavelength, reflectanceSpectra(:, i) .* 100, 'Color', clrs, 'Marker', markers{i}, ...
			'LineWidth', 1.3, 'LineStyle', lineStyle, 'DisplayName', currentCurvName); % plot estimated reflectances
    end
    if (curveN > 1)
        plot(wavelength(peakIdx), reflectanceSpectra(peakIdx, 1) .* 100, 'rx', 'DisplayName', curveNames{1}, 'LineWidth', 1.3); % plot measured reflectance
    end
	hold off
	
	xlabel('Wavelength \lambda (nm)', 'FontSize', 15);
	ylabel('Reflectance %', 'FontSize', 15);
	xlim([400, 700]);
	%figTitle = strjoin({'Comparative plot of Wiener estimation for Sample', simpleSampleName(figTitle)}, ' ');
	%title(figTitle);
    if (curveN > 1)
        legend({curveNames{2:(curveN+1)}, curveNames{1}}, 'Location', 'best', 'FontSize', 15); % 'Orientation','horizontal');
    else
        legend(curveNames{1}, 'Location', 'northwest', 'FontSize', 15); % 'Orientation','horizontal');
    end

    saveOptions.cropBorders = true;
	savePlot(fig, saveOptions);

end

function [sampleName] = simpleSampleName(name)
    attr = strsplit(name);
    sampleName = attr{1};
end