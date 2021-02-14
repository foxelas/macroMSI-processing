function result = applyRowFunc(funcName, varargin)
%APPLYROWFUNC applies a function on each row
%
%   result = applyRowFunc(func, varargin) applies function func with
%   arguments varargin on each row of varargin
%

expectedArgs = nargin(funcName);

rows = size(varargin{1}, 1);
result = zeros(rows, 1);
newVarargin = cell(1, expectedArgs);
for i = 1:rows
    for j = 1:expectedArgs
        newVarargin{j} = varargin{j}(i, :);
    end
    result(i) = funcName(newVarargin{:});
end

end