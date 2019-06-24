function [] = plotSingleIllumination(illumination,fig,saveOptions)
	
	if (nargin < 2)
        illumination = 'd65';
    end
    if (nargin < 2)
        fig = figure;
    else     
        figure(fig);
        clf(fig);
    end
    if (nargin < 3)
        saveOptions.SaveImage = false;
    end
	
	[lambda, illum] = illuminant(illumination);
	plot(lambda, illum/max(illum), 'm', 'LineWidth', 2);
	xlabel('wavelength (nm)')
	xlim([300, 830])
	ylabel('relative spectral power distribution')
	title(strjoin({'The CIE ', upper(illumination) , 'daylight illuminant '}, {' '}))
	savePlot(fig, saveOptions);

end
