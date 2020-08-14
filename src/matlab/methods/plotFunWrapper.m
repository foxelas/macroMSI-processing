function plotFunWrapper(fig, functionName, varargin)

if isnumeric(fig) && ~isempty(fig)
    %disp('Check if no overlaps appear and correct fig is saved.')
    figure(fig);
    clf(fig);
else
    fig = gcf;
end

newVarargin = varargin;
expectedArgs = nargin(functionName);
for i = (length(newVarargin) + 1):(expectedArgs - 1)
    newVarargin{i} = [];
end
newVarargin{length(newVarargin)+1} = fig;
functionName(newVarargin{:});
end