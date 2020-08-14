function recovered = registerImage(moving, tform, newDims)

height = newDims(1);
width = newDims(2);

if ndims(moving) == 2
    outputView = imref2d(newDims);
    recovered = imwarp(moving, tform, 'OutputView', outputView);
elseif ndims(moving) == 3
    channels = size(moving, 1);
    outputView = imref2d([height, width]);
    recovered = zeros([channels, height, width]);
    for i = 1:channels
        recovered(i, :, :) = imwarp(squeeze(moving(i, :, :)), tform, 'OutputView', outputView);
    end

elseif ndims(moving) == 4
    channels = size(moving, 1);
    rgbChannels = size(moving, 4);
    outputView = imref2d([height, width]);
    recovered = zeros([channels, height, width, rgbChannels]);
    for i = 1:channels
        for j = 1:rgbChannels
            recovered(i, :, :, j) = imwarp(squeeze(moving(i, :, :, j)), tform, 'OutputView', outputView);
        end
    end
else
    disp('Unsupported image size.')
end

end