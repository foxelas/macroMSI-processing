function M = getMarkovian(N, rho)

    if (nargin < 2)
        rho =  0.97;
    end
    
    M_row = ones(1,N);
    for i=2:N
        M_row(i) = M_row(i-1) * rho;
    end

    M = zeros(N); %first order Markov process covariance matrix
    select = eye(1, N);
    D = zeros(N,N); % Q x lambdaN
    for i=1:N
        M(i,:) = [ fliplr(M_row(1: i)) , M_row(2:N-i+1)];
        D(i,:) = circshift( select', (i-1)*N); %diagonal with 1 in the position of components
    end
    
end