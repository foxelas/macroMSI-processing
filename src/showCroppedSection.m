function [] = showCroppedSection(idx, x, y, mask, figtitle, savedir, color)

    if ~exist('x', 'var') || ~exist('y', 'var')
        x = idx.Originx;
        y = idx.Originy;
    end

    if (nargin < 4)
        infile = matfile(fullfile('..', '..','input', 'saitama_v3_min_region', 'in.mat'));
        mask = infile.MSIstruct(idx.Index).MaskI;
    end

    if (nargin < 5) || isempty(figtitle)
        figtitle = '';
    end

    if (nargin < 6) || isempty(savedir)
        savedir = '';
        toSave = false;
    else
        toSave = true;
    end

    if (nargin < 7)
        color = 'r';
    end

    datafile = matfile(fullfile('..', '..','input', 'saitama_v3_min_region', 'data.mat'));
    files = {datafile.data(idx.Data).File};
    [~, im] = readMSI(files);
    I = im + 0.3 * cat(3, mask, mask, mask);
    w = warning('off', 'all');

    figure(1);
    hold on
    imshow(I);
    plot(x, y, strcat(color,'*'), 'LineWidth', 2, 'MarkerSize', 5);
    hold off
    %rectangle('Position', [x, y, width, height], 'EdgeColor', color);

    title(figtitle)
    if toSave
        print(fig, savedir, '-djpeg');
    end

    warning(w);

end
