function rmse = rmse(reconstructed, measured)
% Root Mean Square Error
N = size(measured, 2);
rmse = sqrt(((measured - reconstructed) * (measured - reconstructed)')/N);
end