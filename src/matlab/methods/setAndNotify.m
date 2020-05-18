function paramValue = setAndNotify(paramName, paramValue)

    onOffOptions = {'OFF', 'ON'};
    if islogical(paramValue)
        fprintf('Executing with [%s] set to %s.\n',paramName, onOffOptions{paramValue + 1} );
    elseif ischar(paramValue)
        fprintf('Executing with [%s] set to %s.\n',paramName, paramValue );
    elseif isnumeric(paramValue)
        fprintf('Executing with [%s] set to %.2f.\n',paramName, paramValue );
    else 
        warning('Unsupported variable type.\n')
    end 
end

