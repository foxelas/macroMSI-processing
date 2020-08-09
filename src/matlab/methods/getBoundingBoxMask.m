function [bboxMask, bbox] = getBoundingBoxMask(corners, mask)

bbox = getBoundingBox(corners);
[X,Y] = meshgrid(1:size(mask,2), 1:size(mask,1));
bboxMask = mask & X >= corners(1) & X <= corners(2) & Y >= corners(3) & Y <= corners(4) ;

end