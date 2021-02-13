function [x] = getWavelenghts(m, option)
%GETWAVELENGTHS returns the wavelengths 
%   
%   x = getWavelengths(m) returns wavelengths as a vecrtor of wavelengths 
%   x = getWavelengths(m, 'raw') returns wavelengths as a vecrtor of
%   wavelengths
%   x = getWavelengths(m, 'index') returns indexes respective to selected
%   wavelengths 
%

    if nargin < 2 
        option = 'raw';
    end 
    
    switch option 
        case 'raw'
            splitWavelength = getSetting('splitWavelength'); 
            if m == 401 
                x = 380:780;
            elseif m == 36
                x = 380:10:730;
            elseif m == 17
                x = 380:10:splitWavelength;
            elseif m == 19
                x = (splitWavelength + 1):10:730;
            elseif m == 161
                x = 380:splitWavelength;
            elseif m == 240
                x = (splitWavelength + 1):780;
            else 
                error('Unsupported wavelength range');
            end
    
        case 'index'
            x = getWavelenghts(m, 'raw');
            x = x - 380 + 1; 
            
        otherwise 
            error('Unsupported option.')

    end 
    
end 