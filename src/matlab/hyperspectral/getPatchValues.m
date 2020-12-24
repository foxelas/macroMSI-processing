function [actualVals] = getPatchValues(target, colorMasks)

n = size(colorMasks, 3);
v = size(target, 3);
actualVals = zeros(n, v);

if ndims(colorMasks) < 3
    colorMasks = reshape(colorMasks, [size(colorMasks, 1), size(colorMasks, 2), 1]);
end
for k = 1:n
    mask = colorMasks(:, :, k);
    maskedLab = target(any(mask, 2), any(mask, 1), :);
    actualVals(k, :) = mean(reshape(maskedLab, [size(maskedLab, 1) * size(maskedLab, 2), v]));
end

end