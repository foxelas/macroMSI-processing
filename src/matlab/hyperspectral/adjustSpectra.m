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

if (sum(contains(lineNames, 'white 9.5 (.05 D)')) == 0)
    adjustmentMethod = 'noAdjustment';
end 
switch adjustmentMethod
    case 'toRatio'
        white95Idx = strcmp(lineNames, 'white 9.5 (.05 D)');
        white95Val = 0.933;
        alpha = white95Val / max(inSpectra(white95Idx, :));
        fprintf('Values adjusted so that white 9.5 (.05 D) line is assinged to value 0.9 \nwith multiplication by alpha = %.3f \n', alpha);
        outSpectra = inSpectra * alpha;

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
        fprintf('Values adjusted so that white 9.5 (.05 D) line is assinged to value 0.9 \nwith division by alpha = %.3f \n', alpha);
        outSpectra = inSpectra / alpha;
    case 'noAdjustment'
        outSpectra = inSpectra;
        alpha = 1;
    otherwise
        error('Unsupported method');
end

end