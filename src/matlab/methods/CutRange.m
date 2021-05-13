function [newSpectrum, newX] = CutRange(oldSpectrum, x)
%%CUTRANGE removes noisy bands from the spectrum 

ids = x >= 420 & x <= 730;
newX = x(ids);
newSpectrum = oldSpectrum(ids);
end 