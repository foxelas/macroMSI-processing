function [outSpectra, alpha] = adjustSpectra(inSpectra, lineNames, adjustmentMethod)
% ADJUSTSPECTRA multiplies the spectra with a coefficient
%
%     [outSpectra, alpha] = adjustSpectra(inSpectra, lineNames) adjusts
%     spectra values according to 'fixWhiteLevel' adjustmentMethod, so
%     that the white patch is assigned to 90% reflectance
%
%     [outSpectra, alpha] = adjustSpectra(inSpectra, lineNames,
%     adjustmentMethod) adjusts spectra values according to adjustmentMethod
%

if nargin < 3
    adjustmentMethod = 'fixWhiteLevel';
end

switch adjustmentMethod
    case 'fixWhiteLevel'
        white95Idx = strcmp(lineNames, 'white 9.5 (.05 D)');
        white95Val = 0.9;
        
        if size(inSpectra, 2) == 36
            startIdx = 36 - 20;
        elseif size(inSpectra, 2) == 17
            startIdx = 16;
        else
            startIdx = 1;
        end
        alpha = mean(inSpectra(white95Idx, startIdx:end)) / white95Val;
        fprintf('Values adjusted so that white 9.5 (.05 D) line is assinged to value 0.8 \nwith division by alpha = %.3f \n', alpha);
        outSpectra = inSpectra / alpha;
    otherwise
        error('Unsupported method');
end

end