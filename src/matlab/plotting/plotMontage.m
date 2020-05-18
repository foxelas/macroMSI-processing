function [] = plotMontage(left,right,fig) 

    warning('off');
    
    imshowpair(left,right,'montage');
    pause(0.1)
    savePlot(fig);
    
    warning('on');
end

