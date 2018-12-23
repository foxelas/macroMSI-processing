function [coeff, score, latent, explained] = PCA(X)

mu = mean(X);%average of each column
Xcent = bsxfun(@minus, X ,mu);%data centered
[coeff,latent] = svd(cov(Xcent));%Singular value decomposition
latent = diag(latent); 
explained=100*latent/sum(latent);%variances of all individual principal components
score = Xcent * coeff;%each column of score is a Component 
% plot(score(:,1),score(:,2),'r*')%ploting The two principal components
% title(['The two principal components PCA (', num2str(round(sum(explained(1:3)))), '% variance explained)'])
% xlabel('First component')
% ylabel('Second component ')

end

