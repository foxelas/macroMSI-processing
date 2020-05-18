function [] = plotSystemSpectralCharacteristics(plotType, wavelength, sensitivity, illumination,fig)

	%{
		% Arguments  ( wavelength, bandWavelength, illumination )
		% wavelength: 401x1 or 81x1
		% bandWavelength: the names of the bands
		% illumination: the illumination light 401x7
	%}

    
	if contains(lower(plotType), 'sensitivity')
		sensitivityN = size(sensitivity, 2);
		if getSetting('saveInBW')
			color(1,:)=HInt2RGB(1,100); % red, darkest 
			color(2,:)=HInt2RGB(3,64); % green, less dark 
			color(3,:)=HInt2RGB(7,10); % blue cyan, lightest
		else 
			color = hsv(3);
			%color = flip(color);
		end

		n  = max(sensitivity(:));
		
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
		ylabel('Sensitivity [a.u.]', 'FontSize', 15);    
		xlabel('Wavelength (nm)', 'FontSize', 15);
        ylim([0,1]);
        yticks([0, 0.5, 1]);
	end
	
	if (contains(lower(plotType), 'sensitivity') && contains(lower(plotType), 'illumination')); yyaxis right; end

	if contains(lower(plotType), 'illumination')
		if getSetting('saveInBW')
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
        fc = [450, 465, 505, 525, 575, 605, 630];      
		for j = 1:numel(fc)
			plot(wavelength, illumination(:, j), 'DisplayName', [ num2str(fc(j)), ' nm'], 'Color', color(j,:), 'LineWidth', 2);
		end
		plot(wavelength, illumination(:,8), 'DisplayName', 'white', 'LineStyle', ':', 'LineWidth', 3);
		hold off;
		xlabel('Wavelength (nm)', 'FontSize', 15);
		ylabel('Luminous Intensity [mW/sr/m^2]', 'FontSize', 5);
	end
	
	ax = gca;
	ax.FontSize = 15; 
	xlim([400, 700]);  
	xticks([400, 500, 600, 700]);	
	l = legend;
	l.FontSize = 13;
	
    setSetting('plotName', fullfile(getSetting('savedir'), 'general', plotType));
	savePlot(fig);

end
