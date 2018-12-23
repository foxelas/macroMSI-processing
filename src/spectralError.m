function [ error ] = spectralError( r, rhat )
%spectralError computes the spectral rms error of the reflectance
%estimation
    error = ones(1,2);
    
    N = length(r);
    % Root Mean Square Error
    error(1) = sqrt( (r - rhat)' * (r - rhat) / N );
    
    % Normalized Mean Square Error
    error(2) = (r - rhat)' * (r - rhat) / (r' * r); 
    
end

