function [] = plotLBP(lbpImage, figTitle, fig)

%% Plot LBP values on image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im = imagesc(imrotate(lbpImage, 90));
colormap('parula');
if max(lbpImage(:)) > 0.3
    upperlimit = 0.4;
else
    upperlimit = 0.04;
end
im.AlphaData = .8;
c = colorbar;
axis off;
c.Label.FontSize = 15;
c.Label.FontWeight = 'bold';
c.Label.String = 'Texture Descriptor Value';
c.Limits = [0, upperlimit];
c.LimitsMode = 'manual';

%     title(strjoin({figTitle, 'descriptor values'}, ' '),'FontSize', 20);
%     set(gcf, 'Position', get(0, 'Screensize'));

%%add the other textures
setSetting('cropBorders', true);
savePlot(fig);

end
