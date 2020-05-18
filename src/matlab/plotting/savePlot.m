function [] = savePlot( fig, saveOptions)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if (nargin < 1)
    fig = gcf;
end 

if ~exist('saveOptions', 'var')
    saveOptions.saveImages = false;
end
if ~isfield(saveOptions,'saveInHQ')
    saveOptions.saveInHQ = false;
end
if ~isfield(saveOptions,'saveInBW')
    saveOptions.saveInBW = false;
end
if ~isfield(saveOptions,'plotName')
    saveOptions.plotName = [];
end

if (saveOptions.saveImages)
    if (~isempty(saveOptions.plotName))
        filename = strrep(saveOptions.plotName, '.mat', '');

        [filepath,name,~] = fileparts(filename);
        filepathBW = fullfile(filepath, 'bw');    
        mkNewDir(filepath);
        mkNewDir(filepathBW);

        if (saveOptions.saveInHQ)
            filename = fullfile(filepath, strcat(name, '.png'));
            export_fig(filename , '-png','-native', '-nocrop');
            %print(handle, strcat(plotName, '.png'), '-dpng', '-r600');
        else
            filename = fullfile(filepath, name);
            if (isfield(saveOptions, 'cropBorders') && saveOptions.cropBorders)
                export_fig(filename , '-png','-native');
            else
                saveas(fig, filename, 'png');
            end
            if (saveOptions.saveInBW)
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

