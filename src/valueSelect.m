function [G] = valueSelect(A, method)
%% VALUESELECT Selects one appropriate pixel intensity value from
% valueSelect chooses the appropriate values from the image to be used in
% the processing, depending on the 'method'
% method : {'green', 'rms', 'adjusted', 'extended', 'unchanged'}

height = size(A, 2);
width = size(A, 3);

G = zeros(7, height, width);
for row = 1:height
    for col = 1:width
        switch method
            case 'green'
                G(:, row, col) = squeeze(A(:, row, col, 2));
                
            case 'rms'
                G(:, row, col) = squeeze(sqrt(sum(A(:, row, col, :).^2, 4)));
                
            case 'adjusted'
                G(1:2, row, col) = squeeze(A(1:2, row, col, 3)); % blue range
                G(3:5, row, col) = squeeze(A(3:5, row, col, 2)); % green range
                G(6:7, row, col) = squeeze(A(6:7, row, col, 1)); % red range
                
            case 'extended'
                G(1:3, row, col) = squeeze(A(1:3, row, col, 3)); % extended blue range
                G(4:7, row, col) = squeeze(A(3:6, row, col, 2)); % extended green range
                G(8:9, row, col) = squeeze(A(6:7, row, col, 1)); % extended red range
                
            case 'rgb'
                G = A;
                
            case 'nothing'
                G = A;
            otherwise
                error('Unexpected pixelValueSelectionMethod. Could not compute pixel value matrix.')
        end
        
    end
end


end
