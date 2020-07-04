function [] = plotMontage(left, right, figTitle, fig)

hasTitle = true; 
if nargin < 4
    hasTitle = isnumeric(figTitle);
    if ~hasTitle
        fig = figTitle; 
        figTitle = '';
    end 
end 

warning('off');

imshowpair(left, right, 'montage');
if hasTitle 
    title(figTitle);
end 
pause(0.1)
savePlot(fig);

warning('on');
end
