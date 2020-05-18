function [] = plotReconstructionPerformanceBars(performance,lineNames,figTitle,fig,saveOptions)
%% Plot bar performance for reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (nargin < 4)
        fig = figure;
    else        
        figure(fig);
        clf(fig);
    end

    if (nargin < 5)
        saveOptions.SaveImage = false;
    end
    
    if (saveOptions.saveInBW)
        color(1,:)=HInt2RGB(1,100); % red, darkest 
        color(2,:)=HInt2RGB(3,64); % green, less dark 
        color(3,:)=HInt2RGB(7,10); % blue cyan, lightest
    else 
        color = hsv(3);
        %color = flip(color);
    end
    isStacked = true;
    barPlot = plotBars(lineNames, performance * 100, color, figTitle, 'Malignancy', 'NRMSE(%)',...
        {'unfixed', 'fixed', 'sectioned'}, {'Fixing'}, [0, 10], 12, true, isStacked);
    
    savePlot(fig, saveOptions);

end
