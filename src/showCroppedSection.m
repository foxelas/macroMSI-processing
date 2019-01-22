function [] = showCroppedSection(idx, x, y, mask, color)

    w = warning('off', 'all');
    
    if ~exist('x', 'var') || ~exist('y', 'var')
        x = idx.Originx;
        y = idx.Originy;
    end
    
    if (nargin < 5)
        color = 'g*';
    end

    if (nargin < 4)
        infile = matfile(fullfile('..', '..','input', 'saitama_v3_min_region', 'in.mat'));
        msi = infile.MSIStruct(idx.Index,:);
        mask = msi.MaskI;
    end

    dt = matfile(fullfile('..', '..','input', 'saitama_v3_min_region', 'data.mat'));
    datax = dt.data(1,idx.Data);
    files = {datax.File};
    [~, im] = readMSI(files);
    I = im + 0.5 * cat(3, mask, mask, mask);
    plots('cropped', 'Image', I, 'PlotName', idx.Csvid, 'Markers', color, 'Coordinates', [x,y]);

    warning(w);

end
