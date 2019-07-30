function [] = plotMontage(left,right,fig,saveOptions)

    if (nargin < 3)
        fig = figure;
    else 
        figure(fig);
        clf(fig);
    end

    if (nargin < 4)
        saveOptions.SaveImage = false;
    end 

    warning('off');
    
    imshowpair(left,right,'montage');
    saveOptions.cropBorders = true;
    pause(0.1)
    savePlot(fig, saveOptions);
    
    warning('on');
end

