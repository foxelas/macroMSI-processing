function [] = savePlot(fig)
%SAVEPLOT saves the plot shown in figure fig
%   savePlot(2);

saveImages = getSetting('saveImages');

if (saveImages)
    saveInHQ = getSetting('saveInHQ');
    saveInBW = getSetting('saveInBW');
    plotName = getSetting('plotName');
    if (~isempty(plotName))
        filename = strrep(plotName, '.mat', '');

        [filepath, name, ~] = fileparts(filename);
        filepathBW = fullfile(filepath, 'bw');
        mkNewDir(filepath);
        mkNewDir(filepathBW);

        if (saveInHQ)
            filename = fullfile(filepath, strcat(name, '.png'));
            export_fig(filename, '-png', '-native', '-nocrop');
            %print(handle, strcat(plotName, '.png'), '-dpng', '-r600');
        else
            filename = fullfile(filepath, name);
            if getSetting('cropBorders')
                export_fig(filename, '-png', '-native');
            else
                saveas(fig, filename, 'png');
            end
            if (saveInBW)
                filename = fullfile(filepathBW, strcat(name, '.eps'));
                %export_fig(filename , '-eps', '-transparent', '-r900',  '-gray');
            else
                filename = fullfile(filepath, strcat(name, '.eps'));
                %export_fig(filename , '-eps', '-transparent', '-r900',  '-RGB');
            end
        end
    else
        warning('Empty plotname')
    end
end

end
