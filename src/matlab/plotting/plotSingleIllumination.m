function [] = plotSingleIllumination(illumination, fig)

if (nargin < 2)
    illumination = 'd65';
end

[lambda, illum] = illuminant(illumination);
plot(lambda, illum/max(illum), 'm', 'LineWidth', 2);
xlabel('wavelength (nm)')
xlim([300, 830])
ylabel('relative spectral power distribution')
title(strjoin({'The CIE ', upper(illumination), 'daylight illuminant '}, {' '}))
savePlot(fig);

end
