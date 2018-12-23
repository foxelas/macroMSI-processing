function R = markovianMatrix(rho)
%MARKOVIANMATRIX Creates a matrix based on a markovian process with
%correlation rho

    if (nargin < 1)
        rho = 0.985;
    end 

    R_row = ones(1,wavelengthN);
    for i=2:wavelengthN 
        R_row(i) = R_row(i-1) * rho;
    end
    R = zeros(wavelengthN); %first order Markov process covariance matrix
    select = eye(1, wavelengthN);
    D = zeros(wavelengthN,wavelengthN); % Q x lambdaN
    for i=1:wavelengthN
        R(i,:) = [ fliplr(R_row(1: i)) , R_row(2:wavelengthN-i+1)];
        D(i,:) = circshift( select', (i-1)*wavelengthN); %diagonal with 1 in the position of components
    end
    
end

