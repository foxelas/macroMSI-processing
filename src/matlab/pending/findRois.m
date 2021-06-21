contFlag = true;

while contFlag
i=3;
figure(1); imshow(unfixedIaug);
% imshow(fixedIaug); 
% disp('need to select poly');
% v1 = roipoly;
% p = detectMinEigenFeatures(v1); 
% val1 = floor(min(p.Location)); val2 = floor(max(p.Location)); 
% b1 = [val1(1) val2(1) val1(2) val2(2)] 


roiCornerVals = delimread(fullfile(getSetting('systemdir'), getSetting('roiCornerFileName')), ',', {'num', 'text'});
textVals = roiCornerVals.text; 
roiCornerVals = roiCornerVals.num;
[unfixedRoiCorners, ~, ~] = getRoiCorners(unfixedId, roiCornerVals, textVals); 
[fixedRoiCorners, lesion, registrationType] = getRoiCorners(fixedId, roiCornerVals, textVals); 

[unfixedIaug, unfixedBbox] = getCropAreas(unfixedRoiCorners, roiNames, unfixedMask, unfixedWhiteReference);
[fixedIaug, fixedBbox] = getCropAreas(fixedRoiCorners, roiNames, fixedMask, fixedWhiteReference);
figure(2); montage({unfixedIaug, fixedIaug});
v1 =  fixedBbox{i};
figure(3); montage({imcrop(unfixedWhiteReference, unfixedBbox{i}), imcrop(fixedWhiteReference, v1)});
v2 = [v1(1) v1(1)+v1(3)-1 v1(2) v1(2)+v1(4)-1]
end 

%% Functions 
function [roiCorners, lesion, regist] = getRoiCorners(idd, roiVals, txtVals)
idx = find(roiVals(:, 1) == idd);
roiCorners = {roiVals(idx, 2:5), roiVals(idx, 6:9), roiVals(idx, 10:13)};
lesion = txtVals{idx + 1, 1};
regist = strrep(txtVals{idx + 1, 15}, ' ', '');
end 

function [Iaug, bbox] = getCropAreas(roiCorners, roiNames, mask, white)
    roiColors = {'m', 'g', 'c'};
    roiMask = cell(numel(roiNames), 1);
    bbox = cell(numel(roiNames), 1);

    Iaug = white;
    for i = 1:numel(roiNames)
        [roiMask{i}, bbox{i}] = getBoundingBoxMask(roiCorners{i}, mask);
        if sum(isnan(bbox{i})) < 1
            Iaug = insertShape(Iaug, 'rectangle', bbox{i}, 'LineWidth', 5, 'Color', roiColors{i});
        end
    end
end 
