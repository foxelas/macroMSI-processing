function [ coeff,  score, latent, Ipc1, Ipc2] = performPCA( G, names)
%performPCA Performs Principal Component Analysis on a 3D image (RGB or multispectral)
%   data: 3D or 4D image, the first 3 dimensions represent band -
%   length-height. The last band may represent an RGB color value
%   pixelValueSelectionMethod is used to reduce the 4D image to 3D
 
    
    if (ndims(G) == 3 )
        X = reshape(G,size(G,1)*size(G,2), size(G, 3));        
    else 
            X = G;
    end
    
    [coeff, score, latent, ~, explained]= pca(X);
    PCAtransformedData = X*coeff;
    
%     % Calculate eigenvalues and eigenvectors of the covariance matrix
%     covarianceMatrix = cov(X);
%     [V,D] = eig(covarianceMatrix);
%     
%     % Compare ...
%     disp('Compare coeff with eigenvectors')
%     coeff
%     V
%     disp('Compare transformed data with score')
%     score;
%     
%     % The columns of X*coeff are orthogonal to each other. This is shown with ...
%     disp('Correlation between PCs:')
%     corrcoef(PCAtransformedData)
%     
%     % The variances of these vectors are the eigenvalues of the covariance matrix, and are also the output "latent". 
%     disp('Eigenvalues:')
%     var(PCAtransformedData)'
%     latent
%     sort(diag(D),'descend')
    
    fig1 = figure(1);
    plot(latent, '-mx');
%     title(sprintf('#%s: PCA eigenvalues', strrep(sampleName,'_',' ')));
    latent(1:10)
    explained(1:10)

    plots('pca', 2, score, 'PCA Fix', 'lineNames', names)
    plots('pca', 3, score, 'PCA Sample', 'lineNames', names)


    % Visualize both the orthonormal principal component coefficients for each variable and the principal component scores for each observation in a single plot    
% biplot(coeff(:,1:2),'scores',score(:,1:2),'varlables',strtrim(cellstr(num2str(variables'))'));
    %title({sprintf('#%s: Visualize orthonormal PC1&2 for each variable', strrep(sampleName,'_',' ')), '+ the PC scores for each observation'});
    
%     fig2 = figure(2);
%     hold on 
%     %dummy legends
%     
%     %dummy legends
%     idx = arrayfun(@(x) any(strcmp(x, 'unfixed')),{ID.type});
%     scatter(score(idx,1), score(idx,2), 'ro');
%     hold on
%     idx = arrayfun(@(x) any(strcmp(x, 'fixed')),{ID.type}); 
%     scatter(score(idx,1), score(idx,2), 'bx');
%     hold off 
%     title('PCA scores for PC1&2')
%     xlabel('Principal component 1');
%     ylabel('Principal component 2');
%     legend('unfixed', 'fixed');
%     
%     fig3 = figure(3);
%     for i = 1:length(ID)
%         tmp = strsplit(ID(i).id, '_');
%         samples{i} = tmp{2};
%     end
%    samples = unique(samples);
%    
%    cmap = hsv(length(samples));
%     hold on
%     for i=1:length(samples)
%         idx = arrayfun(@(x) (strfind(x, samples{i}) ),{ID.id});
%         idx = cellfun(@(x) ~isempty(x), idx);        
%         scatter(score(idx,1), score(idx,2), [], cmap(i,:));
%     end
%     hold off
%     title('PCA scores for PC1&2')
%     xlabel('Principal component 1');
%     ylabel('Principal component 2');
%     legend(samples)  
%    
%     
%     if ~(isempty(outName))
%         print(fig1, strcat(outName, '_eigenvalues', '.jpg'), '-djpeg');
%         print(fig2, strcat(outName, '_pc1and2_type', '.jpg'), '-djpeg');
%         print(fig3, strcat(outName, '_pc1and2_sample', '.jpg'), '-djpeg');
%     end
    
    if ndims(G) == 3
        Ipc1 = reshape(PCAtransformedData(:,1),size(G,1), size(G,2));
        Ipc2 = reshape(PCAtransformedData(:,2),size(G,1), size(G,2));
    else
        Ipc1 = PCAtransformedData(:,1);
        Ipc2 = PCAtransformedData(:,2);
    end
    
end

