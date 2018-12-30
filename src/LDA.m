% LDA - MATLAB subroutine to perform linear discriminant analysis
% by Will Dwinnell and Deniz Sevis
%
% Use:
% W = LDA(Input,Target,Priors)
%
% W       = discovered linear coefficients (first column is the constants)
% Input   = predictor data (variables in columns, observations in rows)
% Target  = target variable (class labels)
% Priors  = vector of prior probabilities (optional)
%
% Note: discriminant coefficients are stored in W in the order of unique(Target)
%
% Example:
%
% % Generate example data: 2 groups, of 10 and 15, respectively
% X = [randn(10,2); randn(15,2) + 1.5];  Y = [zeros(10,1); ones(15,1)];
%
% % Calculate linear discriminant coefficients
% W = LDA(X,Y);
%
% % Calulcate linear scores for training data
% L = [ones(25,1) X] * W';
%
% % Calculate class probabilities
% P = exp(L) ./ repmat(sum(exp(L),2),[1 2]);
%
%
% Last modified: Dec-11-2010


function [W, scores] = LDA(Input, Target, Priors, Method)

% Determine size of input data
[n, m] = size(Input);
if (n ~= length(Target))
    error('Inconsistent number of observation values and labels');
end

if (nargin < 4)
    Method = '';
end

if strcmp(Method, 'pdf')
    
    class = unique(Target);
    classes = numel(class);
    classIdx = cell(classes, 1);
    
    % class means columns
    mu = zeros(m, classes);
    elements = zeros(classes, 1);
    for i = 1:classes
        classIdx{i} = find(Target == class(i));
        mu(:, i) = mean(Input(classIdx{i}, :));
        elements(i) = numel(classIdx{i});
    end
    globalMu = mean(Input)';
    
    % class variances / class scatter matrices
    S = zeros(classes, m, m);
    for i = 1:classes
        %S(i, :, :) = cov(Input(classIdx{i},:));
        S(i, :, :) = (Input(classIdx{i}, :)' - mu(:, i)) * (Input(classIdx{i}, :)' - mu(:, i))';
    end
    
    % within-class scatter matrix
    Sw = squeeze(sum(S)); % sw = squeeze(s(1, :, :)) + squeeze(s(2, :,:));
    
    % between-class scatter matrix
    Sb = (mu - globalMu) * diag(elements) * (mu - globalMu)';
    
    % eigenvalues
    invSwSb = Sw \ Sb;
    [V, ~] = svd(invSwSb);
    
    % the projection vector
    W = V(:, :); %eigenvector of the maximum eigenvalue % W2 = inv(Sw) * (mu(:,1) - mu(:,2));
    scores = Input * W;
    
    %         transMean = W' * mu;
    
elseif strcmp(Method, 'kardi')
    
    class = unique(Target);
    classes = numel(class);
    classIdx = cell(classes, 1);
    
    % class means columns
    mu = zeros(m, classes);
    elements = zeros(classes, 1);
    for i = 1:classes
        classIdx{i} = find(Target == class(i));
        mu(:, i) = mean(Input(classIdx{i}, :));
        elements(i) = numel(classIdx{i});
    end
    globalMu = mean(Input);
    
    centeredX = Input - globalMu;
    
    % class variances / pooled variance matrices
    S = zeros(classes, m, m);
    poolCov = zeros(m, m);
    
    if ~isempty(Priors)
        % Use the user-supplied priors
        PriorProb = Priors;
    else
        % Use the sample probabilities
        PriorProb = elements / n;
    end
    
    for i = 1:classes
        S(i, :, :) = centeredX(classIdx{i}, :)' * centeredX(classIdx{i}, :) / numel(classIdx{i});
        poolCov = poolCov + squeeze(S(i, :, :)*numel(classIdx{i})/n);
    end
    invPoolCov = inv(poolCov);
    
    W = zeros(classes, m+1);
    for i = 1:classes
        temp = mu(:, i)' / invPoolCov;
        W(i, :) = [-0.5 * temp * mu(:, i) + log(PriorProb(i)), temp];
    end
    
    scores = [ones(n, 1), Input] * W';
    
else
    % Discover and count unique class labels
    ClassLabel = unique(Target);
    k = length(ClassLabel);
    
    % Initialize
    nGroup = NaN(k, 1); % Group counts
    GroupMean = NaN(k, m); % Group sample means
    PooledCov = zeros(m, m); % Pooled covariance
    W = NaN(k, m+1); % model coefficients
    
    % Loop over classes to perform intermediate calculations
    for i = 1:k
        % Establish location and size of each class
        Group = (Target == ClassLabel(i));
        nGroup(i) = sum(double(Group));
        
        % Calculate group mean vectors
        GroupMean(i, :) = mean(Input(Group, :));
        
        % Accumulate pooled covariance information
        PooledCov = PooledCov + ((nGroup(i) - 1) / (n - k)) .* cov(Input(Group, :));
    end
    
    % Assign prior probabilities
    if ~isempty(Priors)
        % Use the user-supplied priors
        PriorProb = Priors;
    else
        % Use the sample probabilities
        PriorProb = nGroup / n;
    end
    
    % Loop over classes to calculate linear discriminant coefficients
    for i = 1:k
        % Intermediate calculation for efficiency
        % This replaces:  GroupMean(g,:) * inv(PooledCov)
        Temp = GroupMean(i, :) / PooledCov;
        
        % Constant
        W(i, 1) = -0.5 * Temp * GroupMean(i, :)' + log(PriorProb(i));
        
        % Linear
        W(i, 2:end) = Temp;
    end
    
    scores = [ones(n, 1), Input] * W';
    
    % Housekeeping
    clear Temp
end
end


% EOF
