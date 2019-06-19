function [] = plotMontage(left,right,fig,saveOptions)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    if (nargin < 3)
        fig = figure;
    else 
        figure(fig);
    end

    if (nargin < 4)
        saveOptions.SaveImage = false;
    end 


    imshowpair(left,right,'montage');
    
    savePlot(fig, saveOptions);
end

