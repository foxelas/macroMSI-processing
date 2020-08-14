function poiId = groupId2poiId(groupId, ID)
[~, idxs] = find([ID.Group] == groupId, 1);
poiId = ID(idxs).Index;
end