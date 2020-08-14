function [] = savePlot(fig)
%SAVEPLOT saves the plot shown in figure fig
%   savePlot(2);

saveImages = getSetting('saveImages');

if (saveImages)
    saveInHQ = getSetting('saveInHQ');
    saveInBW = getSetting('saveInBW');
    plotName = getSetting('plotName');
    cropBorders = getSetting('cropBorders');
    saveEps = getSetting('saveEps');

    if (~isempty(plotName))
        filename = strrep(plotName, '.mat', '');

        [filepath, name, ~] = fileparts(filename);
        filepathBW = fullfile(filepath, 'bw');
        mkNewDir(filepath);
        mkNewDir(filepathBW);

        filename = fullfile(filepath, strcat(name, '.png'));
        if (cropBorders)
            export_fig(filename, '-png', '-native', '-transparent');
        else
            if (saveInHQ)
                export_fig(filename, '-png', '-native', '-nocrop');
                %print(handle, strcat(plotName, '.png'), '-dpng', '-r600');
            else
                saveas(fig, filename, 'png');
            end
        end
        if (saveEps)
            namext = strcat(name, '.eps');
            if (saveInBW)
                filename = fullfile(filepathBW, namext);
                export_fig(filename, '-eps', '-transparent', '-r900', '-gray');
            else
                filename = fullfile(filepath, namext);
                export_fig(filename, '-eps', '-transparent', '-r900', '-RGB');
            end
        end
    else
        warning('Empty plotname')
    end
end

end
