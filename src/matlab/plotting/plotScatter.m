function plotScatter(xx, yy, maskAll, maskRois, maskColors, xlab, ylab)

xx = xx(maskAll);
yy = yy(maskAll);
for i = 1:size(maskRois, 1)
    m = maskRois{i};
    maskRois{i} = logical(m(maskAll));
end

figure;
hold on;
scatter(xx, yy);
for i = 1:size(maskRois, 1)
    scatter(xx(maskRois{i}), yy(maskRois{i}), [], maskColors{i});
end 
xlabel(xlab);
ylabel(ylab);

hold off;
end 