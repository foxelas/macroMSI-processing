function [] = plotMap(I, type, fig)

c = imagesc(I);
c.Parent.Visible = 'off';
colorbar;
switch type
    case 'opticalDensityMelanin'
        title('OD630nm - Melanin Map');
    case 'opticalDensityHemoglobin'
        title('OD575nm - 1.15 OD630nm - Hemoglobin Map');
    case 'deepMelanin'
        title('Deep Melanin');
    case 'totalMelanin'
        title('Total Melanin');
end

setSetting('cropBorders', true);
setSetting('plotName', fullfile(getSetting('savedir'), getSetting('map'), ...
    type, getSetting('outName')));

savePlot(fig);

end
