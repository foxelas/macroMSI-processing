function [G] = raw2msi(A, method)
%% RAW2MSI Selects one appropriate pixel intensity value from
% raw2msi chooses the appropriate values from the image to be used in
% the processing, depending on the 'method'
% method : {'green', 'rms', 'adjusted', 'extended', 'unchanged', 'max'}


[bands, height, width, ~] = size(A);

switch method
	case 'green'
		G = squeeze(A(:, :, :, 2));
		
	case 'rms'
		G = squeeze(sqrt(sum(A.^2, 4)));
    
    case 'max'
        G = squeeze(max(A.^2, [], 4));
		
	case 'adjusted'
		G = zeros(bands, height, width);
		G(1:2, :, :) = squeeze(A(1:2, :, :, 3)); % blue range
		G(3:5, :, :) = squeeze(A(3:5, :, :, 2)); % green range
		G(6:7, :, :) = squeeze(A(6:7, :, :, 1)); % red range
		
	case 'extended'
		bands = 9;
		G = zeros(bands, height, width);
		G(1:3, :, :) = squeeze(A(1:3, :, :, 3)); % extended blue range
		G(4:7, :, :) = squeeze(A(3:6, :, :, 2)); % extended green range
		G(8:9, :, :) = squeeze(A(6:7, :, :, 1)); % extended red range
		
	case 'rgb'
		G = A;
		
	case 'nothing'
		G = A;
	otherwise
		error('Unexpected pixelValueSelectionMethod. Could not compute pixel value matrix.')
end
        


end
