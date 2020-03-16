k = 20;
msiType = 'extended'; %'extended'; % 'max';

% for all pixels without bg removal 
removebg = false; 
[msi, whiteReference, specimenMask, height, width, channels] = getImage(k, options, msiType, removebg);
columns = reshape(msi, channels, width * height)'; %30k pixels x 7 variable
[coeff,latent,explained]  = doPixelPCA(columns);

% for all pixels with bg removal 
removebg = true; 
[msi, whiteReference, specimenMask, height, width, channels] = getImage(k, options, msiType, removebg);
columns = reshape(msi, channels, width * height)'; %30k pixels x 7 variable
[coeff,latent,explained]  = doPixelPCA(columns);

% for only the pixels of the speciment 
fgColumn = reshape(specimenMask, 1, width*height);
columns = reshape(msi, channels, width * height)'; %30k pixels x 7 variable
columns = columns(fgColumn, :);
[coeff,latent,explained]  = doPixelPCA(columns);


function [coeff,latent,explained]  = doPixelPCA(columns)

[coeff,score,latent,tsquared,explained,mu] = pca(columns,'NumComponents',2, 'Centered', true);

coeff
latent
explained

end 