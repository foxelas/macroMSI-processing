function [fileConditions] = getFileConditions(content, target, id)
%%GETFILECONDITIONS returns the conditions necessary for finding the
%%filename of the file to be read 
%
%   Usage: 
%   fileConditions = getFileConditions(content, target)

    fileConditions = {getSetting('configuration'), [], ...
        getSetting('integrationTime'), target, getSetting('dataDate')};
    
    if nargin > 1 
        fileConditions = {getSetting('configuration'), content, ...
        getSetting('integrationTime'), target, getSetting('dataDate')};
    end 
    
    if nargin > 2 
        fileConditions = [fileConditions, id];
    end 
end 