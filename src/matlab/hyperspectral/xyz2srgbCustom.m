% ================================================
% *** FUNCTION xyz2srgb
% ***
% *** function [RGB] = xyz2srgb(XYZ)
% *** computes 8-bit sRGB from XYZ
% *** XYZ is n by 3 and in the range 0?1
% *** see also srgb2xyz
function [RGB] = xyz2srgbCustom(XYZ)
if (size(XYZ, 2) ~= 3)
    disp('XYZ must be n by 3');
    return;
end
M = [3.2404542, -1.5371385, -0.4985314; ...
    -0.9692660, 1.8760108, 0.0415560; ...
    0.0556434, -0.2040259, 1.0572252];
RGB = (M * XYZ')';
RGB(RGB < 0) = 0;
RGB(RGB > 1) = 1;
index = (RGB <= 0.00304);
RGB = RGB + (index) .* (12.92 * RGB);
RGB = RGB + (1 - index) .* (1.055 * RGB.^(1 / 2.4) - 0.055);

end