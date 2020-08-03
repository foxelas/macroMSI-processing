function [regMsi, regWhiteReference, regSpecimenMask] = registerAllRelated(msi, whiteReference, specimenMask, tform, newDims)
    regMsi = registerImage(msi, tform, newDims);
    regWhiteReference = permute(registerImage(permute(whiteReference, [3, 1, 2]), tform, newDims), [2, 3, 1]);
    regSpecimenMask = registerImage(specimenMask, tform, newDims);
end 