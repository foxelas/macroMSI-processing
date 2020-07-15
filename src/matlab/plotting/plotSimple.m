function [] = plotSimple(I, figTitle, fig)

hasTitle = true; 
if nargin < 3
    hasTitle = isnumeric(figTitle);
    if ~hasTitle
        fig = figTitle; 
        figTitle = '';
    end 
end 

warning('off');

imshow(I, []);
if hasTitle 
    title(figTitle);
end 
pause(0.1)
savePlot(fig);

warning('on');
end