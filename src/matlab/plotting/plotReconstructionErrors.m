function [] = plotReconstructionErrors(errors, figTitle,fig,saveOptions)
	if (nargin < 3)
        fig = figure;
    else     
        figure(fig);
        clf(fig);
    end
    if (nargin < 4)
        saveOptions.SaveImage = false;
    end
	
	if contains(figTitle, 'nmse', 'IgnoreCase', true)
		avge = [errors.avgrmse];
		maxe = [errors.maxrmse];
		mine = [errors.minrmse];
		stde = [errors.stdrmse];
		labely = 'NMSE';
	else
		avge = [errors.avgnmse];
		maxe = [errors.maxnmse];
		mine = [errors.minnmse];
		stde = [errors.stdnmse];
		labely = 'RMSE';
	end
	
	if contains(figTitle, 'System') 
		pvsms = unique({errors.pixelValueSelectionMethod}, 'stable');
	end
	
	if contains(figTitle, 'Noise') 
		pvsms = unique({errors.noiseType}, 'stable');
	end
	
	if contains(figTitle, 'Matrix')
		smms = unique({errors.smoothingMatrixMethod}, 'stable');
	end 
	labelx = 'Wiener estimation variation';

	
	marker = ['o', 's', 'd', '^', '*', 'h', 'p', 'v', '<', '+', '>','o', 's', 'd', '^', '*'];
	if ((length(smms) > 1) && (length(pvsms) > 1))
		pvsmColor = parula(length(pvsms)+1);
		pvsmColor = pvsmColor(1:length(pvsms),:);
	else
		pvsmColor = jet(length(pvsms));
	end
%         pvsmColor = pvsmColor(2:9,:); % remove yellow and black color 
	
	hold on
	% Generate dummy info for plot handles "h"
	h = zeros(length(pvsms), 1);
	for i = 1:length(pvsms)
		h(i) = plot([NaN, NaN], 'Color', pvsmColor(i, :), 'Marker', marker(i), 'DisplayName', pvsms{i});
	end
	if contains(figTitle, 'avg', 'IgnoreCase', true)
		i = i + 1;
		h(i) = plot([NaN, NaN], 'ko', 'DisplayName', 'Mean and Std');
	end
	if contains(figTitle, 'max', 'IgnoreCase', true)
		i = i + 1;
		h(i) = plot([NaN, NaN], 'k-.', 'DisplayName', 'Max value');
	end
	if contains(figTitle, 'min', 'IgnoreCase', true)
		i = i + 1;
		h(i) = plot([NaN, NaN], 'k:', 'DisplayName', 'Min value');
	end
	% end dummy data for legend entries
	
	for i = 1:length(pvsms)
		for j = 1:length(smms)
			if contains(figTitle, 'avg', 'IgnoreCase', true)
				errorbar(j, avge((i - 1)*length(smms)+j), stde((i - 1)*length(smms)+j) ./ size(stde,2) , ...
					'Color', pvsmColor(i, :), 'LineWidth', 1, 'Marker', marker(i), 'MarkerSize', 10);
			end
		end
		if contains(figTitle, 'max', 'IgnoreCase', true)
			plot(1:length(smms), maxe, 'Color', pvsmColor(i, :), 'LineStyle', '-.');
		end
		if contains(figTitle, 'min', 'IgnoreCase', true)
			plot(1:length(smms), mine, 'Color', pvsmColor(i, :), 'LineStyle', ':');
		end
	end
	hold off
	
	grid on
	grid minor
	xlabel(labelx);
	xlim([0, numel(smms) + 1]);
	xticks(1:numel(smms));
	xticklabels(strrep(smms, '_', ' '));
	xtickangle(45);
	ylabel(labely);
	ax = gca;
	ax.YRuler.Exponent = 0;
	legend(h, 'Location', 'best');
	title('Comparison of Mean and Standard Error Values for various Estimation Configurations')
	set(gcf, 'Position', get(0, 'Screensize'));
		
	saveOptions.saveInHQ = true;
		
    savePlot(fig, saveOptions);

end