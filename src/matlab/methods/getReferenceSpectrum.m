function [reference] = getReferenceSpectrum(n)
    if nargin < 1
        n = 81;
    end 
    if n == 81
        load(fullfile('parameters', 'referenceSpectrum.mat'), 'referenceSpectrum81');
        reference = referenceSpectrum81';
    else 
        load(fullfile('parameters', 'referenceSpectrum.mat'), 'referenceSpectrum401');
        reference = referenceSpectrum401';
    end 
end 