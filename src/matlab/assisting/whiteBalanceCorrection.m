close all; 
locationDir = '..\..\..\..\..\..\mspi\saitama_v2\white\unfixed\right';
whiteFilename = 'P2270082.tif';
msiFilenames = {'P2270083.tif', 'P2270084.tif', 'P2270085.tif', 'P2270086.tif', 'P2270087.tif', 'P2270088.tif', 'P2270089.tif'};
allFilenames = [whiteFilename, msiFilenames];

%% Test various color correction techniques on white light image 
A = imread(fullfile(locationDir, whiteFilename));
A_sRGB = lin2rgb(A);
warning('off','images:initSize:adjustingMag')
figure(1);montage({A,A_sRGB})
title('Original, minimally processed image before and after gamma correction')

%     % Use roipoly to create a mask from a polygon drawn manually. 
%     % Click to add vertices, then right-click and select "Create Mask" to return.
%     title('Draw a polygon around the chart') 
%     mask_chart = roipoly; 
%     mask_chart = imdilate(mask_chart,ones(7));
load('saved parameters\mask_chart.mat', 'mask_chart');

% Use the provided estimate of ROI center coordinates.
x = [1262 1778 2274 2746 3210 3690];
y = [2498 2490 2458 2450 2450 2434];

x = round(x);
y = round(y);
r = mean(diff(x)) / 2 * 0.60;
r = floor(r);

mask = false(size(A,1), size(A,2));
for k = 1:6
    mask(y(k)-r:y(k)+r,x(k)-r:x(k)+r) = true;
end

mask_eroded = imerode(mask, strel('disk',5));

mask_clipped = (A == intmax(class(A))) | (A == intmin(class(A)));
mask_clipped = mask_clipped(:,:,1) | mask_clipped(:,:,2) | mask_clipped(:,:,3);

mask_patches = mask_eroded & ~mask_clipped;

A_patches = imoverlay(A_sRGB,mask_patches);

figure(2);imshow(A_patches)
title('The selected pixels are highlighted in yellow')

patches_R = A(:,:,1);
patches_G = A(:,:,2);
patches_B = A(:,:,3);
patches_R = patches_R(mask_patches);
patches_G = patches_G(mask_patches);
patches_B = patches_B(mask_patches);

patches_R = im2double(patches_R);
patches_G = im2double(patches_G);
patches_B = im2double(patches_B);

illuminant_groundtruth = [mean(patches_R) mean(patches_G) mean(patches_B)];

illuminant = [0.066 0.1262 0.0691];
figure(3);
plot3([0 1],[0 1],[0,1],'LineStyle',':','Color','k')
hold on
plot3(...
    [0 illuminant_groundtruth(1)/norm(illuminant_groundtruth)], ... % Red
    [0 illuminant_groundtruth(2)/norm(illuminant_groundtruth)], ... % Green
    [0 illuminant_groundtruth(3)/norm(illuminant_groundtruth)], ... % Blue
    'Marker','.', 'MarkerSize',10)
hold on
plot3( ...
    [0 illuminant(1)/norm(illuminant)], ... % Red
    [0 illuminant(2)/norm(illuminant)], ... % Green
    [0 illuminant(3)/norm(illuminant)], ... % Blue
    'Marker','.', 'MarkerSize',10)
xlabel('R')
ylabel('G')
zlabel('B')
title('Illuminants in RGB space')
xlim([0 1])
ylim([0 1])
zlim([0 1])
view(28, 36)
legend('achromatic line', 'ground truth illuminant', 'estimated illuminant')
grid on
axis equal

illuminant_gw1 = illumgray(A, 1.5, 'Mask', ~mask_chart);
err_gw1 = colorangle(illuminant_gw1, illuminant_groundtruth);
disp(['Angular error for Gray World with percentiles=[0 0]: ' num2str(err_gw1)])
B_gw1 = chromadapt(A, illuminant_gw1, 'ColorSpace', 'linear-rgb');

B_gw1_sRGB = lin2rgb(B_gw1);
figure(4);
imshow(B_gw1_sRGB)
title('White balanced image using Gray World with percentiles=[0 0]')

illuminant_ch1 = illumpca(A, 5, 'Mask', ~mask_chart);
err_ch1 = colorangle(illuminant_ch1, illuminant_groundtruth);
disp(['Angular error for Cheng with percentage=5: ' num2str(err_ch1)])
B_ch1 = chromadapt(A, illuminant_ch1, 'ColorSpace', 'linear-rgb');
B_ch1_sRGB = lin2rgb(B_ch1);
figure(5);
imshow(B_ch1_sRGB)
title('White balanced image using Cheng with percentage=5')

illuminant_ch2 = illumpca(A, 'Mask', ~mask_chart);
err_ch2 = colorangle(illuminant_ch2, illuminant_groundtruth);
disp(['Angular error for Cheng with percentage=3.5: ' num2str(err_ch2)])
B_ch2 = chromadapt(A, illuminant_ch2, 'ColorSpace', 'linear-rgb');
B_ch2_sRGB = lin2rgb(B_ch2);
figure(6);
imshow(B_ch2_sRGB)
title('White balanced image using Cheng with percentile=3.5')

illuminant_wp1 = illumwhite(A, 0, 'Mask', ~mask_chart);
err_wp1 = colorangle(illuminant_wp1, illuminant_groundtruth);
disp(['Angular error for White Patch with percentile=0: ' num2str(err_wp1)])
B_wp1 = chromadapt(A, illuminant_wp1, 'ColorSpace', 'linear-rgb');
B_wp1_sRGB = lin2rgb(B_wp1);

figure(7);
imshow(B_wp1_sRGB)
title('White balanced image using White Patch Retinex with percentile=0')


illuminant_wp2 = illumwhite(A, 1, 'Mask', ~mask_chart);
err_wp2 = colorangle(illuminant_wp2, illuminant_groundtruth);
disp(['Angular error for White Patch with percentile=1: ' num2str(err_wp2)])
B_wp2 = chromadapt(A, illuminant_wp2, 'ColorSpace', 'linear-rgb');
B_wp2_sRGB = lin2rgb(B_wp2);

figure(8);
imshow(B_wp2_sRGB)
title('White balanced image using White Patch Retinex with percentile=1')

save('saved parameters\color_correction.mat', 'illuminant_gw1')