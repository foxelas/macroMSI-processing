function [error] = finderror(xx, datafile, measured, id, sensitivity, illumination, options)
    width = 5;
    height = 5;

    x = xx(1);
    y = xx(2);

    g = readMSI( datafile, fix(x), fix(y), width, height);
    [~, error, ~] = reflectanceEstimation( g, measured, id, options);

end
