mask = unfixedMask;
 
cornersHb = [316, 382, 242, 295];   %bboxHb = [316 242 382-316 295-242];
cornersNorm = [159, 252, 167, 214];    %bboxNorm = [159 154 252-159 225-154];
cornersMel = [398, 440, 83, 137];   %bboxMel = [398 83 440-398 137-83];
[maskHb, bboxHb] = getBoundingBoxMask(cornersHb, mask);
[maskNorm, bboxNorm] = getBoundingBoxMask(cornersNorm, mask);
[maskMel, bboxMel] = getBoundingBoxMask(cornersMel, mask);


Iaug = insertShape(I,'rectangle', bboxMel,'LineWidth',5, 'Color', 'c');
Iaug = insertShape(Iaug,'rectangle', bboxHb,'LineWidth',5, 'Color', 'm');
Iaug = insertShape(Iaug,'rectangle', bboxNorm,'LineWidth',5, 'Color', 'g');
imshow(Iaug)


setSetting('saveImages', false);

metrics = getImageSimilarity(imcrop(hbMap1, bboxHb), imcrop(hbMap2, bboxHb), [], unfixedId, fixedId, 'unfixed', 'fixed');

