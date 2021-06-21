function [hasGpu] = hasGPU()
%%HASGPU informs whether there is gpu available 

    v = dbstack;
    if numel(v) > 1
        parentName = v(2).name;
    else
        parentName = 'none';
    end
    isFirst = contains(parentName, 'initialization');
    if isFirst
        %pcName = char(java.net.InetAddress.getLocalHost.getHostName);
        %if stcmp(pcName, 'GPU-PC2') == 0
        if length(ver('parallel')) == 1
            setSetting('pcHasGPU', true);
        end 
    end 
    hasGpu = getSetting('pcHasGPU');
end 