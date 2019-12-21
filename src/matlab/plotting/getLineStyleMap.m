function [lineStyleMap] = getLineStyleMap(style)

    if (nargin  < 1 )
        style = 'class';
    end 
    
    switch style
        case 'class'
            key = {'Benign', 'Atypical', 'Malignant'};
            value = {'go', 'ms', 'r^'}; 
        case 'type'
            key = {'Unfixed', 'Fixed', 'Sectioned'};
            value = {'go', 'ms', 'r^'}; 
        case 'sample'
            key =   {'0037', '0045', '0053', '0059', '0067', '9913', '9933', '9940', '9949', '9956'};
            value =	{'o', 'x', 'd', '^', '*', 'h', 'p', 'v', 's', '<', '+', '>'};
    end
            
    lineStyleMap = containers.Map(key, value);
end

