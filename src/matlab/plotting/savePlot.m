function [] = savePlot( fig )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if (nargin < 1)
    fig = gcf;
else     
    figure(fig);
    clf(fig);
end

saveImages = getSetting('saveImages');
saveInHQ = getSetting('saveInHQ');
saveInBW = getSetting('saveInBW');
plotName = getSetting('plotName');

if (saveImages)
    if (~isempty(plotName))
        filename = strrep(plotName, '.mat', '');

        [filepath,name,~] = fileparts(filename);
        filepathBW = fullfile(filepath, 'bw');    
        mkNewDir(filepath);
        mkNewDir(filepathBW);

        if (saveInHQ)
            filename = fullfile(filepath, strcat(name, '.png'));
            export_fig(filename , '-png','-native', '-nocrop');
            %print(handle, strcat(plotName, '.png'), '-dpng', '-r600');
        else
            filename = fullfile(filepath, name);
            if getSetting('cropBorders')
                export_fig(filename , '-png','-native');
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

