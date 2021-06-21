function plotScatterHbMel(hpMap, melMap, mask, roiBboxes)

[m, n] = size(mask);
maskCol = reshape(mask, [m*n, 1]);
hbMapCol = reshape(hpMap, [m*n,1]);
melMapCol = reshape(melMap, [m*n, 1]);

colors = {'m', 'g', 'c'};
colorMasks = cell(numel(roiBboxes), 1);
for i = 1:numel(roiBboxes)
    roiBbox = roiBboxes{i};
    xStart = roiBbox(1);
    yStart = roiBbox(2);
    xEnd = xStart + roiBbox(3) - 1;
    yEnd = yStart + roiBbox(4) - 1;
    colorMask = zeros(m, n);
    colorMask(yStart:yEnd, xStart:xEnd) = 1; 
    colorMasks{i} = reshape(colorMask, [m*n, 1]);
    
end

plotScatter(melMapCol, hbMapCol, maskCol, colorMasks, colors, 'Total Melanin', 'Total Hemoglobin');

end 