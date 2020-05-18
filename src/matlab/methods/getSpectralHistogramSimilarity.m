function [val] = getSpectralHistogramSimilarity(h1, h2)
    val = 0; 
    for i = 1:size(h1, 1)
        val = val + sum((h1 - h2).^2 ./ (h1 + h2));
    end
    mn = sum(h1(1,:));
    val = 1 / mn * sum(val);
end 
