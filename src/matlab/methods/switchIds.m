function [ID] = switchIds(ID, data, arr)
%     SWITHCIDS switches the ids depending on the raw image illumination settings 
% 
%     Usage:
%     ID = switchIds(ID, data, 'brights')
%     ID = switchIds(ID, data, 'darks')

groups = findgroups([ID.MsiID]);
for g = 1:max(groups)
    gIdxs = find(groups == g);
    [ID(gIdxs).MsiID] = deal(arr(g).MsiID);
    [ID(gIdxs).Path] = deal(arr(g).Path);
    [ID(gIdxs).Directory] = deal(arr(g).Directory);
    [ID(gIdxs).Filename] = deal(arr(g).Filename);
    [ID(gIdxs).RgbID] = deal(find([data.MsiID] == arr(g).MsiID, 1));
end
end