function [x] = getWavelengths(m, option)
%GETWAVELENGTHS returns the wavelengths
%
%   x = getWavelengths(m) returns wavelengths as a vector of wavelengths
%   x = getWavelengths(m, 'raw') returns wavelengths as a vecrtor of
%   wavelengths
%   x = getWavelengths(m, 'index') returns indexes respective to selected
%   wavelengths
%   x = getWavelengths(m, 'babel') returns indexes respective to selected
%   wavelengths for babel standard spectra
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
        x = getWavelengths(m, 'raw');
        x = x - 380 + 1;
        
    case 'babel'
        if m == 36
            x = 1:36;
        elseif m == 17
            x = 1:17;
        elseif m == 19
            x = 18:36;
        else
            error('Unsupported wavelengths range');
        end
        
    otherwise
        error('Unsupported option.')
        
end

end