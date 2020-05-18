function [normMsi] = getNormMSI(msi)
%% getNormMSI returns normalizations of each channel by removing its mean

    if ndims(msi) == 3 
        normMsi = zeros(size(msi));
        for i = 1:size(msi,1)
            normMsi(i, :, :) = squeeze(msi(i, :, :)) - mean(mean(squeeze(msi(i, :, :))));
        end 
    elseif ndims(msi) == 2
        normMsi = zeros(size(msi));
        for i = 1:size(msi,2)
            normMsi(:,i) = squeeze(msi(:, i)) - mean(squeeze(msi(:, i)));
        end 
    else
        disp('Not supported dimensions')
    end 
end