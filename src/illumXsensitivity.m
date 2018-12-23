function [ H ] = illumXsensitivity( E, S, pixelValueSelectionMethod )
%illumXsensitivity Computes the system matrix H using illumination and
%sensitivity
%   Input args
%   E: illumination , spectrumdim x MSdim e.g. 401x7
%   S: camera sensitivity, spectrumdim x filterdim, e.g. 401x3, RGB
%   pixelValueSelectionMethod: either of {'green', 'rms', 'adjusted', 'extended', 'unchanged'}
%   Output args
%   H: system matrix, MSdim x spectrumdim e.g. 7x401

    if (size(S, 2) > 3)
        if strcmp( pixelValueSelectionMethod, 'extended' )
            pixelValueSelectionMethod = 'fullbandext';
        else 
            pixelValueSelectionMethod = 'fullband';
        end
    end
    
    if (size(E,1) ~= size(S,1))
        
    end
    
    switch pixelValueSelectionMethod
        case 'green'
            H = E' * diag(S(:,2));

        case  'rms'
            H = E' * diag( sqrt(sum(S.^2,2)) );

        case 'adjusted'
            H(1:2,:) = E(:,1:2)' * diag(S(:,3)); %  blue range
            H(3:5,:) = E(:, 3:5)' * diag(S(:,2)); %  green range
            H(6:7,:) = E(:, 6:7)' * diag(S(:,1)); %  red range

        case 'extended'
            H(1:3,:) = E(:,1:3)' * diag(S(:,3)); % extended blue range
            H(4:7,:) = E(:, 3:6)' * diag(S(:,2)); % extended green range
            H(8:9,:) = E(:, 6:7)' * diag(S(:,1)); % extended red range

        case 'unchanged' %estimation using RGB images     
            H = zeros(3, length(S));
            lightBands = [6,7; 3,5; 1,2]; %RGB light bands 
            for k=1:3
                Eadj = rms( E( :, lightBands(k,1):lightBands(k,2)), 2 )'; % average the light at certain wavelength bands 
                H(k, 1:length(S)) = Eadj .* S(:,k)';
            end
         
        case 'rgb'
            H(1,:) = E' * diag(S(:,3)); % extended blue range
            H(2,:) = E' * diag(S(:,2)); % extended green range
            H(3,:) = E' * diag(S(:,1)); % extended red range
            
        case 'fullband'
            H = (E .* S)';
        
        case 'fullbandext'
            E = [ E(:, 1:3), E(:, 3:6), E(:, 6:7)];
            S = [ S(:, 1:3), S(:, 3:6), S(:, 6:7)];
            H = (E .* S)'; 
            
        otherwise 
            error('Unexpected pixelValueSelectionMethod. Could not compute system matrix.')
    end 

end

