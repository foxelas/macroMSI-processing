function [difMsi, avgMsi] = getDifMSI(msi, type)
%% Get subimages of differences among band subimages 
    if strcmp(type, 'toAverage')
        channels = size(msi, 1);
        avgMsi =  mean(msi,1);
        avgMsiRep = repmat(avgMsi, channels, 1, 1);   
        difMsi = bsxfun(@(x,y) x - y, msi, avgMsiRep); 
    elseif strcmp(type, 'toNextBand')
        [channels, m, n] = size(msi);
        difMsi = zeros(channels - 1, m, n);
        for i = 1:(channels - 1)
            difMsi(i, :, :) = squeeze(msi(i, :, :) - msi(i + 1, :, :));
        end 
        avgMsi = [];
    else 
        disp('Not supported method.');
    end    
end 