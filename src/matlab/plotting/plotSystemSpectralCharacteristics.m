function [] = plotSystemSpectralCharacteristics(plotType, wavelength, sensitivity, illumination,fig,saveOptions)

	%{
		% Arguments  ( wavelength, bandWavelength, illumination )
		% wavelength: 401x1 or 81x1
		% bandWavelength: the names of the bands
		% illumination: the illumination light 401x7
	%}
	
	if (nargin < 5)
        fig = figure;
    else     
        figure(fig);
        clf(fig);
    end
    if (nargin < 6)
        saveOptions.SaveImage = false;
    end
    
	
	if contains(lower(plotType), 'sensitivity')
		sensitivityN = size(sensitivity, 2);
		if (saveInBW)
			color(1,:)=HInt2RGB(1,100); % red, darkest 
			color(2,:)=HInt2RGB(3,64); % green, less dark 
			color(3,:)=HInt2RGB(7,10); % blue cyan, lightest
		else 
			color = hsv(3);
			%color = flip(color);
		end

		if strcmp(plotType, 'normSensitivity'); n  = max(sensitivity(:)); else; n = 1; end
		
		hold on;
		if (sensitivityN == 7)
			fnames = cellstr(num2str(defaultfc'));
		else
			fnames = {'red'; 'green'; 'blue'};
		end
		for j = 1:sensitivityN
			plot(wavelength, sensitivity(:, j) / n, 'Color', color(j,:), 'LineWidth', 3, 'DisplayName', fnames{j});
		end
		hold off;
		if strcmp(plotType, 'normSensitivity'); ylabel('Normalized Sensitivity', 'FontSize', 17); else; ylabel('Sensitivity', 'FontSize', 17); end     
		xlabel('Wavelength \lambda (nm)', 'FontSize', 15);
	end
	
	if (contains(lower(plotType), 'sensitivity') && contains(lower(plotType), 'illumination')); yyaxis right; end

	if contains(lower(plotType), 'illumination')
		if (saveInBW)
			color(7,:)=HInt2RGB(1,100); % red, darkest 
			color(6,:)=HInt2RGB(7,82); % cyan, less dark 
			color(5,:)=HInt2RGB(3,64); % green, less dark 
			color(4,:)=HInt2RGB(5,53); %
			color(3,:)=HInt2RGB(9,46); % magenta, less dark 
			color(2,:)=HInt2RGB(2,28); % orange, less dark 
			color(1,:)=HInt2RGB(7,10); % blue cyan, lightest
		else 
			color = jet(7);
		end

		hold on;
		for j = 1:numel(fc)
			plot(wavelength, illumination(:, j) * 10, 'DisplayName', [ num2str(fc(j)), ' nm'], 'Color', color(j,:), 'LineWidth', 2);
		end
		plot(wavelength, illumination(:,8) * 10, 'DisplayName', 'white', 'LineStyle', ':', 'LineWidth', 3);
		hold off;
		xlabel('Wavelength \lambda (nm)', 'FontSize', 15);
		ylabel('Luminous Intensity', 'FontSize', 17);
	end
	
	if (contains(lower(plotType), 'sensitivity') && contains(lower(plotType), 'illumination'))
		title('Overlap of illumination and camera sensitivity spectrum')
	elseif contains(lower(plotType), 'normsensitivity') 
		title('Normalized Camera Sensitivities', 'FontSize', 15);
	elseif contains(lower(plotType), 'sensitivity') 
		title('Camera Sensitivities', 'FontSize', 15);
	elseif contains(lower(plotType), 'illumination') 
		title('Normalized Luminous Intensity', 'FontSize', 15);
	end

	ax = gca;
	ax.FontSize = 15; 
	xlim([400, 700]);  
	xticks([400, 500, 600, 700]);
	yticks([0, 0.5, 1]);
	l = legend;
	l.FontSize = 13;
	
	savePlot(fig, saveOptions);

end
