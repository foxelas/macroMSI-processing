function [outSpectrum] = pads(inSpectrum, options)
%PADS adds or removes padding to spectrum array 
% 
%   [outSpectrum] = pads(inSpectrum) adds padding to inSpectrum 
%   [outSpectrum] = pads(inSpectrum, 'add') adds padding to inSpectrum
%   [outSpectrum] = pads(inSpectrum, 'del') removes padding from inSpectrum
%

    if nargin < 2 
        options = 'add';
    end 
    
    m = length(inSpectrum);
    
    switch options 
        case 'add'
            x = getWavelengths(m, 'index');
            if m > 100 
                outSpectrum = zeros(401, 1);
            else 
                outSpectrum = zeros(36, 1);
            end 
            outSpectrum(x) = inSpectrum;
            
        case 'del'
            if m == 36 
                cutoff = 15; 
            else 
                cutoff = 150;
            end 

            if isempty(nonzeros(inSpectrum(1:cutoff)))
                idStart = find(inSpectrum, 1);
                outSpectrum = inSpectrum(idStart:end);
                
            elseif isempty(nonzeros(inSpectrum((end-cutoff):end)))
                idEnd = find(inSpectrum, 1, 'last');
                outSpectrum = inSpectrum(1:idEnd);
                
            else 
                outSpectrum = inSpectrum;
            end
            
        otherwise 
            error('Unsupported options');
    end 
end 