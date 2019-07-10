function [W, score, latent, explained] = dimensionReduction(method, input, labels, priors, targetDimension)
%% DIMENSIONREDUCTION provides data projections for dimension reduction
% Available methods
% 'PCA' Principal Component Analysis, MATLAB built-in
% 'PCA b' Variation of PCA
% 'LDA' Linear Dimension Analysis
% 'LDA b' Variation of LDA
% 'PCALDA' Combination or PCA and LDA 
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
% [W, score, latent, explained] = dimensionReduction('PCA', X);
% [W, score] = dimensionReduction('PCA', X, [], [], 2);
% [W, score] = dimensionReduction('LDA', X, labels);
% [W, score] = dimensionReduction('LDA', X, labels, [0.3, 0.7], 2 );

[n, m] = size(input); %n observations, m variables
    
if (nargin < 3)
    labels = [];
end

if (nargin < 4)
    priors = [];
end

if (nargin < 5)
    targetDimension = m;
end

latent = zeros(10,1);
explained = zeros(10,1);
    
if contains(method, 'LDA')
    % Determine size of input data
    if (n ~= length(labels))
        error('Inconsistent number of observation values and labels');
    end
    
    % Available classes
    classes = unique(labels);
    C = numel(classes);

    if iscell(labels)
        textLabels = labels;
        textClasses = classes;
        labels = zeros(size(labels));
        classes = 1:C;
        for i = 1:C
            labels(strcmp(textLabels, textClasses{i})) = i;
        end
    end
    % Class element counts
    nGroup = zeros(C, 1); 
    for i = 1:C
        nGroup(i) = sum(double(labels == classes(i)));
    end

    if isempty(priors)
        priors = nGroup / n;
    end

    if (numel(priors)~= C)
        error('Inconsistent number of classes and prior probabilities.')
    end
end


switch method
    case 'PCA b'
        [W, score, latent, ~, explained] = pca(input, 'Centered', true);
    
    case 'PCA'
        mu = mean(input);
        Xcent = bsxfun(@minus, input ,mu);
        [W,latent] = svd(cov(Xcent));
        latent = diag(latent); 
        explained=100*latent/sum(latent);
        score = Xcent * W;
        
    case 'LDA'
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
        [W, latent] = svd(invSwSb); %eigenvector of the maximum eigenvalue 
        score = input * W;
        
    case 'LDA b'
        mu = zeros(C, m); % Group sample means
        poolCov = zeros(m, m); % Pooled covariance
        W = zeros(C, m+1);

        for i = 1:C
            group = (labels == classes(i));
            
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
        score = [ones(n, 1), input] * W';

    case 'PCALDA'
        [~, score, ~, ~] = dimensionReduction('PCA', input, labels);
        [W, score] = dimensionReduction('LDA', score(:, 1:10), labels);

    otherwise
        error('Not implemented dimension reduction method')
end

targetDimension = min(targetDimension, size(W,2));
W = W(:,1:targetDimension);

end

