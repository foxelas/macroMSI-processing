function [lineColorMap] = getLineColorMap(style)

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
end

lineColorMap = containers.Map(key, value);
end
