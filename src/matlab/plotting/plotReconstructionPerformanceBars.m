function [] = plotReconstructionPerformanceBars(performance, lineNames, figTitle, fig)

%% Plot bar performance for reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if getSetting('saveInBW')
    color(1, :) = HInt2RGB(1, 100); % red, darkest
    color(2, :) = HInt2RGB(3, 64); % green, less dark
    color(3, :) = HInt2RGB(7, 10); % blue cyan, lightest
else
    color = hsv(3);
    %color = flip(color);
end
isStacked = true;
plots(fig, @plotBars, lineNames, performance*100, color, figTitle, 'Malignancy', 'NRMSE(%)', ...
    {'unfixed', 'fixed', 'sectioned'}, {'Fixing'}, [0, 10], 12, true, isStacked);

savePlot(fig);

end
