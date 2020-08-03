function [] = plotScaledImage(I, figTitle, limits, fig)

hasTitle = true; 
if nargin < 3 || length(limits) < 2 
    limits = [];
end
if nargin < 4
    hasTitle = isnumeric(figTitle);
    if ~hasTitle
        fig = figTitle; 
        figTitle = '';
    end 
end 

warning('off');

imshow(I, limits);
if hasTitle 
    title(figTitle);
end 
pause(0.1)
savePlot(fig);

warning('on');
end