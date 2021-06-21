function [bbox] = getBoundingBoxMask(corners)
%     GETBOUNDINGBOX returns the settings for a bounding box from the mask indexes
%
%     Usage
%     bbox = getBoundingBox(corners)
%
%     Example
%     corners = [316, 382, 242, 295];
%     bbox = getBoundingBox(corners);
%     returns bboxHb = [316 242 382-316 295-242];

bbox = [corners(1), corners(3), corners(2) - corners(1) + 1, corners(4) - corners(3) + 1];

end