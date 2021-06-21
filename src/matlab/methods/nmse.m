function nmse = nmse(reconstructed, measured)
% Normalized Mean Square Error
nmse = (measured - reconstructed) * (measured - reconstructed)' / (measured * reconstructed');
end