function paramValue = setAndNotify(paramName, paramValue)

    onOffOptions = {'OFF', 'ON'};
    if islogical(paramValue)
        fprintf('Executing with [%s] set to %s',paramName, onOffOptions{paramValue + 1} );
    elseif ischar(paramValue)
        fprintf('Executing with [%s] set to %s',paramName, paramValue );
    elseif isnumeric(paramValue)
        fprintf('Executing with [%s] set to %.2f',paramName, paramValue );
    else 
        warning('Unsupported variable type.')
    end 
end

