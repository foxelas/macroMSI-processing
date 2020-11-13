function result = applyFuncOnRows(actualVals, expectedVals, func)

for i = 1:size(expectedVals, 1)
    result(i) = func(actualVals(i, :), expectedVals(i, :));
end

end