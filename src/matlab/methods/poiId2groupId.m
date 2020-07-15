function groupId = poiId2groupId(poiId, ID)
    [~, idxs] = find([ID.Index] == poiId, 1);
    groupId = ID(idxs).Group;
end