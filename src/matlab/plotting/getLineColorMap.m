function [lineColorMap] = getLineColorMap(style, names)
%     GETLINECOLORMAP returns a linecolor map based on the style
%
%     Usage:
%     [lineColorMap] = getLineColorMap('class')

if (nargin < 1)
    style = 'class';
end

switch style
    case 'class'
        key = {'Benign', 'Atypical', 'Malignant'};
        value = {'g', 'm', 'r'};
    case 'type'
        key = {'Unfixed', 'Fixed', 'Sectioned'};
        value = {'g', 'm', 'r'};
    case 'sample'
        key = {'0037', '0045', '0053', '0059', '0067', '9913', '9933', '9940', '9949', '9956'};
        value = jet(10);
    case 'custom'
        key = names;
        value = jet(length(names));
    otherwise 
        key = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'};
        value = jet(10);
end

if size(value,1) == 1
    value = {value};
end
lineColorMap = containers.Map(key, value);
end
