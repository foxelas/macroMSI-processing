function [W, score, latent, explained] = dimensionReduction(method, input, varargin)
%% DIMENSIONREDUCTION provides data projections for dimension reduction
% Available methods
% 'pca' Principal Component Analysis
% 'lda' Linear Dimension Analysis
% 'pcalda' Combination or PCA and LDA 
% 
% Use:
% [W, score, latent, explained] = dimensionReduction(method, input, varargin)
% 
% W         = discovered projections
% score     = projection scores
% latent    = (only for PCA)
% explained = (only for PCA)
% 
% Example: 
% 
% [W, score, latent, explained] = dimensionReduction('pca', X);
% [W, score] = dimensionReduction('pca', X, 'TargetDimension', 2);
% [W, score] = dimensionReduction('lda', X, labels);
% [W, score] = dimensionReduction('lda', X, labels, 'TargetDimension', 2, 'Priors', [0.3, 0.7]);
    
p = inputParser;
p.addRequired(p, 'method', @ischar);
p.addRequired(p, 'input', @isnumeric);
p.addOptional(p, 'labels', [] );
p.addParameter(p, 'Priors', [], @(x) isnumeric(x));
p.addParameter(p, 'TargetDimension', [], @(x) isnumeric(x));
p.parse(p, method, input, varargin{:})

method = p.Results.method;
input = p.Results.input;
labels = p.Results.labels;
priors = p.Results.Priors;
targetDimension = p.Results.TargetDimension;

if contains(method, 'lda')
    % Determine size of input data
    [n, ~] = size(input);
    if (n ~= length(labels))
        error('Inconsistent number of observation values and labels');
    end

    % Available classes
    classes = unique(labels);
    C = numel(classes);

    % Class element counts
    nGroup = zeros(C, 1); 
    for i = 1:C
        group = (labels == classes(i));
        nGroup(i) = sum(double(group));
    end

    if isempty(priors)
        priors = nGroup / n;
    end

    if (numel(priors)~= C)
        error('Inconsistent number of classes and prior probabilities.')
    end

    latent = [];
    explained = [];
end


switch method
    case 'pca'
        [W, score, latent, explained] = PCA(input, 'Centered', true);

    case 'pca b'
        mu = mean(X);
        Xcent = bsxfun(@minus, X ,mu);
        [coeff,latent] = svd(cov(Xcent));
        latent = diag(latent); 
        explained=100*latent/sum(latent);
        score = Xcent * coeff;

    case 'lda'
        classIdx = cell(C, 1);

        % class means columns
        mu = zeros(m, C);
        elements = zeros(C, 1);
        for i = 1:C
            classIdx{i} = find(labels == classes(i));
            mu(:, i) = mean(input(classIdx{i}, :));
            elements(i) = numel(classIdx{i});
        end
        globalMu = mean(input)';

        % class variances / class scatter matrices
        S = zeros(C, m, m);
        for i = 1:C
            %S(i, :, :) = cov(Input(classIdx{i},:));
            S(i, :, :) = (input(classIdx{i}, :)' - mu(:, i)) * (input(classIdx{i}, :)' - mu(:, i))';
        end

        % within-class scatter matrix
        Sw = squeeze(sum(S)); % sw = squeeze(s(1, :, :)) + squeeze(s(2, :,:));

        % between-class scatter matrix
        Sb = (mu - globalMu) * diag(elements) * (mu - globalMu)';

        % eigenvalues
        invSwSb = Sw \ Sb;
        [V, latent] = svd(invSwSb);

        % the projection vector
        W = V(:, :); %eigenvector of the maximum eigenvalue 
        % W2 = inv(Sw) * (mu(:,1) - mu(:,2));
        score = Input * W;
        %         transMean = W' * mu;
        
    case 'lda b'
        mu = zeros(C, m); % Group sample means
        poolCov = zeros(m, m); % Pooled covariance
        W = zeros(C, m+1);

        for i = 1:C
            % Calculate group mean vectors
            mu(i, :) = mean(input(group, :));

            % Accumulate pooled covariance information
            poolCov = poolCov + ((nGroup(i) - 1) / (n - C)) .* cov(input(group, :));
        end

        for i = 1:C
            tmp = mu(i, :) / poolCov;
            % Constant
            W(i, 1) = -0.5 * tmp * mu(i, :)' + log(priors(i));
            % Linear
            W(i, 2:end) = tmp;
        end
        score = [ones(n, 1), Input] * W';

    case 'pcalda'
        [~, score, ~, ~] = dimensionReduction('pca', input, labels);
        [W, score] = dimensionReduction('lda', score(:, 1:10), labels);

    otherwise
        error('Not implemented dimension reduction method')
end

if ~isempty(targetDimension)
    W = W(:,1:n);
    score = score(:,1:n);
end

end

